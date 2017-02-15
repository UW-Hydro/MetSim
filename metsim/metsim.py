

"""
Handles the synchronization of multiple processes for MetSim
"""

import os
import struct
import time as tm
from getpass import getuser
from multiprocessing import Pool

import numpy as np
import pandas as pd
import xarray as xr
# from netCDF4 import Dataset

from metsim.methods import mtclim
import metsim.constants as cnst

references = '''Thornton, P.E., and S.W. Running, 1999. An improved algorithm for estimating incident daily solar radiation from measurements of temperature, humidity, and precipitation. Agricultural and Forest Meteorology, 93:211-228.
Kimball, J.S., S.W. Running, and R. Nemani, 1997. An improved method for estimating surface humidity from daily minimum temperature. Agricultural and Forest Meteorology, 85:87-98.
Bohn, T. J., B. Livneh, J. W. Oyler, S. W. Running, B. Nijssen, and D. P. Lettenmaier, 2013a: Global evaluation of MT-CLIM and related algorithms for forcing of ecological and hydrological models, Agr. Forest. Meteorol., 176, 38-49, doi:10.1016/j.agrformet.2013.03.003.'''

now = tm.ctime(tm.time())
user = getuser()

attrs = {'pet': {'units': 'mm d-1', 'long_name': 'potential evaporation',
                 'standard_name': 'water_potential_evaporation_flux'},
         'prec': {'units': 'mm d-1', 'long_name': 'precipitation',
                  'standard_name': 'precipitation_flux'},
         'swrad': {'units': 'W m-2', 'long_name': 'shortwave radiation',
                   'standard_name': 'surface_downwelling_shortwave_flux'},
         'lwrad': {'units': 'W m-2', 'long_name': 'longwave radiation',
                   'standard_name': 'surface_downwelling_longwave_flux'},
         't_max': {'units': 'C', 'long_name': 'maximum daily air temperature',
                   'standard_name': 'daily_maximum_air_temperature'},
         't_min': {'units': 'C', 'long_name': 'minimum daily air temperature',
                   'standard_name': 'daily_minimum_air_temperature'},
         '_global': {'conventions': '1.6', 'title': 'Output from MetSim',
                     'institution': 'University of Washington',
                     'source': 'metsim.py',
                     'history': 'Created: {0} by {1}'.format(now, user),
                     'references': references,
                     'comment': 'no comment at this time'}}

class MetSim(object):
    """
    MetSim handles the distribution of jobs that write to a common file
    by launching muliple processes and queueing up their writeback so that
    work can be done while IO is happening.
    """

    # Class variables
    process_handles = []
    methods = {'mtclim': mtclim}
    params = {
        "method": '',
        "domain": '',
        "out_dir": '',
        "start": '',
        "stop": '',
        "time_step": '',
        "out_format": '',
        "in_format": None,
        "base_elev": 0,
        "t_max_lr": '',
        "t_min_lr": '',
        "out_vars": ['temp', 'prec', 'shortwave', 'longwave',
                     'wind', 'vapor_pressure', 'rel_humid']
    }

    def __init__(self, params: dict):
        """
        Constructor
        """
        # Record parameters
        self.update(params)
        MetSim.params['dates'] = pd.date_range(params['start'], params['stop'])
        self.method = MetSim.methods[params['method']]
        self.output = None
        self.met_data = None
        self.ready = False

    def load(self, job_list, domain):
        """Preprocess the input data"""
        # Get the necessary information from the domain
        self.domain_variables = ['mask', 'elev', 'lat', 'lon']
        self.domain = xr.open_dataset(domain).rename(
                MetSim.params['domain_vars'])
        self.lat = self.domain['lat']
        self.lon = self.domain['lon']
        self.mask = self.domain['mask']
        self.elev = self.domain['elev']
        self.domain_shape = self.mask.shape
        self.domain_dims = self.mask.dims

        # Input preprocessing 
        in_preprocess = {"ascii" : self.ascii_in_preprocess,
                         "binary" : self.binary_in_preprocess,
                         "netcdf" : self.netcdf_in_preprocess}
        self.met_data = in_preprocess[MetSim.params['in_format']](job_list)
        
        # Output preprocessing 
        func = {"ascii": self.ascii_out_preprocess,
                "netcdf": self.netcdf_out_preprocess}
        func[MetSim.params['out_format']]()

        self.ready = True

    def launch(self, job_list):
        """Farm out the jobs to separate processes"""
        # Load in the data
        self.load(job_list, self.params['domain'])
        nprocs = MetSim.params['nprocs']
        self.pool = Pool(processes=nprocs)

        iinds, jinds = np.nonzero(self.mask.values)
        # TODO: fix the hard coded 10
        locations = np.array_split(list(zip(iinds, jinds)), nprocs * 10)

        # Do the forcing generation and disaggregation if required
        disagg = int(MetSim.params['time_step']) < cnst.MIN_PER_DAY
        status = []
        for loc_chunk in locations:
            stat = self.pool.apply_async(
                wrap_run,
                args=(self.method.run, loc_chunk, self.met_data, disagg),
                callback=self.unpack_results)
            status.append(stat)
        self.pool.close()

        # Check that everything worked
        [stat.get() for stat in status]
        self.pool.join()
        self.write()

    def unpack_results(self, results):
        """Put results into the master dataset"""
        for result in results:
            locs, df = result
            i = locs['lat']
            j = locs['lon']
            for varname in MetSim.params['out_vars']:
                self.output[varname].values[:, i, j] = df[varname].values

    def run(self, locations):
        """
        Kicks off the disaggregation and queues up data for IO
        """
        results = []
        disagg = int(MetSim.params['time_step']) < cnst.MIN_PER_DAY

        for i, j in locations:
            print("Processing {} {}".format(i, j))
            locs = dict(lat=i, lon=j)
            ds = self.met_data.isel(x=j, y=i, **locs)  # fix this indexing
            lat = ds['lat'].values
            elev = ds['elev'].values
            df = ds.drop(['lat', 'lon', 'elev']).to_dataframe()
            df = MetSim.method.run(df, MetSim.params, elev=elev, lat=lat, disagg=disagg)
            results.append((locs, df))

        self.unpack_results(results)

    def find_elevation(self, lat: float, lon: float) -> float:
        """ Use the domain file to get the elevation """
        i_idx = np.abs(self.lat - lat).argmin()
        j_idx = np.abs(self.lon - lon).argmin()
        return self.domain['elev'][i_idx, j_idx]

    def update(self, new_params):
        """Updates the global parameters dictionary"""
        MetSim.params.update(new_params)

    def netcdf_in_preprocess(self, filename):
        """Get the extent and load the data"""
        # Get the information for splitting up the job
        in_forcing = xr.open_dataset(filename).rename(MetSim.params['in_vars'])

        ilist, jlist = np.nonzero(self.mask.values)

        self.i_idx = ilist
        self.j_idx = jlist

        self.locations = list(zip(ilist, jlist))

        print('found %d locations' % len(self.locations))

        if (in_forcing['lat'].ndim == 1) and (in_forcing['lon'].ndim == 1):
            self.lats = in_forcing['lat'].values[ilist]
            self.lons = in_forcing['lon'].values[jlist]
        elif (in_forcing['lat'].ndim == 2) and (in_forcing['lon'].ndim == 2):
            self.lats = in_forcing['lat'].values[ilist, jlist]
            self.lons = in_forcing['lon'].values[ilist, jlist]
        else:
            raise ValueError('lat/lon dimenssions do not match')

        in_forcing.close()
        # Now read in the file properly
        met_data = self.read(filename)
        return met_data

    def ascii_in_preprocess(self, job_list):
        """Process all files to find spatial extent"""
        # Binary forcing files have naming format $NAME_$LAT_$LON
        sets = [os.path.basename(fpath).split("_") for fpath in job_list]
        self.locations = [(float(s[1]), float(s[2])) for s in sets]
        self.lats = np.unique(sorted([float(s[1]) for s in sets]))
        self.lons = np.unique(sorted([float(s[2]) for s in sets]))
        self.lat_idx = {str(lat): i for i, lat in enumerate(self.lats)}
        self.lon_idx = {str(lon): j for j, lon in enumerate(self.lons)}

        # Creates the master dataset which will be used to parallelize
        met_data = xr.Dataset(coords={'time' : MetSim.params['dates'],
                                      'lon' : self.lons,
                                      'lat' : self.lats},
                              attrs={'n_days' : len(MetSim.params['dates'])})
        shape = (len(MetSim.params['dates']), len(self.lats), len(self.lons))

        # Create the empty variables
        met_data['elev'] = (('lat', 'lon'),
                           np.full((len(self.lats), len(self.lons)), np.nan))
        for var in MetSim.params['in_vars']:
            met_data[var] = (('time', 'lat', 'lon'),np.full(shape, np.nan))

        # Fill in the data
        for job in job_list:
            _, lat, lon = os.path.basename(job).split("_")
            ds = self.read(job)
            met_data['elev'][self.lat_idx[lat], self.lon_idx[lon]] = (
                    self.find_elevation(float(lat), float(lon)))
            for var in MetSim.params['in_vars']:
                met_data[var][:, self.lat_idx[lat], self.lon_idx[lon]] = ds[var]
        return met_data

    def binary_in_preprocess(self, job_list):
        """Process all files to find spatial extent"""
        # Binary forcing files have naming format $NAME_$LAT_$LON
        sets = [os.path.basename(fpath).split("_") for fpath in job_list]
        self.locations = [(float(s[1]), float(s[2])) for s in sets]
        self.lats = np.unique(sorted([float(s[1]) for s in sets]))
        self.lons = np.unique(sorted([float(s[2]) for s in sets]))

        # Creates the master dataset which will be used to parallelize
        met_data = xr.Dataset(coords={'time': MetSim.params['dates'],
                                      'lon': self.lons,
                                      'lat': self.lats},
                              attrs={'n_days': len(MetSim.params['dates'])})
        shape = (len(MetSim.params['dates']), len(self.lats), len(self.lons))

        # Create the empty variables
        met_data['elev'] = xr.Variable(
            (('lat', 'lon'),
             np.full((len(self.lats), len(self.lons)), np.nan)))
        for var in MetSim.params['in_vars'].values():
            met_data[var] = xr.Variable(('time', 'lat', 'lon'),
                                        np.full(shape, np.nan))

        # Fill in the data
        for job in job_list:
            _, lat, lon = os.path.basename(job).split("_")
            ds = self.read(job)
            met_data['elev'][self.i_idx[lat], self.j_idx[lon]] = (
                    self.find_elevation(float(lat), float(lon)))
            for var in MetSim.params['in_vars'].values():
                met_data[var][:, self.i_idx[lat], self.j_idx[lon]] = ds[var]
        return met_data

    def ascii_out_preprocess(self):
        """Dummy function"""
        pass

    def netcdf_out_preprocess(self):
        """Initialize the output file"""
        print("Initializing netcdf...")
        out_dir = MetSim.params['out_dir']
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        self.output = self.domain.copy(deep=True)
        self.output.attrs = attrs['_global']
        self.output_filename = os.path.join(out_dir, "forcing.nc")

        # Number of timesteps
        disagg = int(MetSim.params['time_step']) < cnst.MIN_PER_DAY
        if disagg:
            delta = pd.Timedelta('1 days')
        else:
            delta = pd.Timedelta('0 days')
        times = pd.date_range(MetSim.params['start'],
                              MetSim.params['stop'] + delta,
                              freq="{}T".format(MetSim.params['time_step']))
        n_ts = len(times)

        shape = (n_ts, ) + self.domain_shape
        dims = ('time', ) + self.domain_dims
        coords = {'time': times, **self.domain.mask.coords}
        for varname in MetSim.params['out_vars']:
            self.output[varname] = xr.DataArray(
                data=np.full(shape, np.nan),
                coords=coords, dims=dims,
                name=varname, attrs=attrs.get(varname, {}),
                encoding=None)

    def write(self):
        """
        Dispatch to the right function based on the configuration given
        """
        dispatch = {
                'netcdf': self.write_netcdf,
                'ascii': self.write_ascii
                }
        dispatch[MetSim.params.get('out_format', 'netcdf').lower()]()

    def write_netcdf(self):
        """Write out as NetCDF to the output file"""
        print("Writing netcdf...")
        self.output.to_netcdf(self.output_filename,
                              unlimited_dims=['time'],
                              encoding={'time': {'dtype': 'f8'}})

    def write_ascii(self, data: dict):
        """Write out as ASCII to the output file"""
        print("Writing ascii...")
        out_dir = MetSim.params.get('out_dir', './results/')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        for l in self.output.keys():
            lat, lon = l.split("_")
            out_file = '_'.join(["forcing", str(lat), str(lon)])
            data[l].to_csv(os.path.join(out_dir, out_file), sep='\t')

    def read(self, fpath: str) -> xr.Dataset:
        """
        Dispatch to the right function based on the file extension
        """
        dispatch = {
                'binary'   : self.read_binary,
                'netcdf'    : self.read_netcdf,
                'ascii'   : self.read_ascii
                }
        return dispatch[MetSim.params['in_format']](fpath)

    def read_binary(self, fpath: str) -> xr.Dataset:
        """ Reads a binary forcing file (VIC 4 format) """
        dates = pd.date_range(MetSim.params['start'], MetSim.params['stop'])
        n_days = len(dates)
        precip = []  # Short unsigned int
        t_max = []  # Short int
        t_min = []  # Short int
        wind = []  # Short int

        # Pack these for nicer syntax in the loop
        var_name = [precip, t_max, t_min, wind]
        scale = [40.0, 100.0, 100.0, 100.0]

        # Data types referred to: 'H' - unsigned short ; 'h' - short
        types = ['H', 'h', 'h', 'h']
        with open(fpath, 'rb') as f:
            i = 0
            points_read = 0
            points_needed = 4*n_days
            while points_read != points_needed:
                bytes = f.read(2)
                if bytes:
                    # Get correct variable and data type with i,
                    # then unpack & scale
                    var_name[i].append(
                        struct.unpack(types[i], bytes)[0] / scale[i])
                    i = (i + 1) % 4
                    points_read += 1
                else:
                    break

        # Binary forcing files have naming format $NAME_$LAT_$LON
        param_list = os.path.basename(fpath).split("_")
        params = {"name": param_list[0],
                  "lat": float(param_list[1]),
                  "lon": float(param_list[2]),
                  "n_days": int(n_days)}
        MetSim.params.update(params)
        params['elev'] = [[self.find_elevation(params['lat'], params['lon'])]]
        df = xr.Dataset({"prec": (['time'], precip),
                         "t_min": (['time'], t_min),
                         "t_max": (['time'], t_max),
                         "wind": (['time'], wind),
                         "elev": (['lon', 'lat'], params['elev']),
                         "day_of_year": (['time'], dates.dayofyear)},
                        coords={'lon': [params['lon']],
                                'lat': [params['lat']],
                                'time': dates},
                        attrs={'n_days': params['n_days']})
        return df

    def read_ascii(self, fpath:str) -> xr.Dataset:
        """Read in an ascii forcing file"""
        dates = pd.date_range(MetSim.params['start'], MetSim.params['stop'])
        ds = pd.read_table(fpath, header=None, names=MetSim.params['in_vars'])
        ds.index = dates
        ds = xr.Dataset.from_dataframe(ds)
        return ds


    def read_netcdf(self, fpath:str) -> xr.Dataset:
        """
        Read in a NetCDF file and add elevation information
        """
        dates = pd.date_range(MetSim.params['start'], MetSim.params['stop'])
        n_days = len(dates)
        ds = xr.open_dataset(fpath)
        ds = ds.sel(time=slice(MetSim.params['start'], MetSim.params['stop']))
        ds.rename(MetSim.params['in_vars'], inplace=True)

        # Add elevation and day of year data
        ds['elev'] = xr.Variable(['y', 'x'], self.elev.values)
        ds['day_of_year'] = xr.Variable(['time'], dates.dayofyear)

        # Update the configuration
        MetSim.params.update({"n_days": n_days})
        return ds


def wrap_run(func, loc_chunk, met_data, disagg):
    # this is wrapped so we can return a tuple of locs and df
    results = []
    for i, j in loc_chunk:
        print("Processing {} {}".format(i, j))
        locs = dict(lat=i, lon=j)
        ds = met_data.isel(x=j, y=i, **locs)  # fix this indexing
        lat = ds['lat'].values
        elev = ds['elev'].values
        df = ds.drop(['lat', 'lon', 'elev']).to_dataframe()
        df = func(df, MetSim.params, elev=elev, lat=lat, disagg=disagg)
        results.append((locs, df))
    return results

