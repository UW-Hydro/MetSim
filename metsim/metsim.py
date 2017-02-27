

"""
Handles the synchronization of multiple processes for MetSim
"""

import os
import struct
import itertools
import time as tm
from getpass import getuser
from multiprocessing import Pool

import numpy as np
import pandas as pd
import xarray as xr

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
        "sw_prec_thresh": 0.0,
        "mtclim_swe_corr": False,
        "lw_cloud": 'cloud_deardorff',
        "lw_type": 'prata',
        "tdew_tol": 1e-3,
        "tmax_daylength_fraction": 0.67,
        "out_vars": ['temp', 'prec', 'shortwave', 'longwave',
                     'wind', 'vapor_pressure', 'rel_humid']
    }

    def __init__(self, params: dict):
        """
        Constructor
        """
        # Record parameters
        MetSim.params.update(params)
        MetSim.params['dates'] = pd.date_range(params['start'], params['stop'])

        self.output = None
        self.met_data = None
        self.ready = False

    def load(self, job_list):
        """Load the necessary datasets into memory"""
        # Get the necessary information from the domain
        domain = self.params['domain']
        self.domain = xr.open_dataset(domain).rename(
                MetSim.params['domain_vars'])
        self.lat = self.domain['lat']
        self.lon = self.domain['lon']
        self.mask = self.domain['mask']
        self.elev = self.domain['elev']
        self.domain_shape = self.mask.shape
        self.domain_dims = self.mask.dims
        self.disagg = int(MetSim.params['time_step']) < cnst.MIN_PER_DAY
        self.method = MetSim.methods[self.params['method']]
        ilist, jlist = np.nonzero(np.nan_to_num(self.mask.values))

        self.i_idx = ilist
        self.j_idx = jlist

        self.locations = list(zip(ilist, jlist))
        print('found %d locations' % len(self.locations))

        out_dir = MetSim.params['out_dir']
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        self.output = self.domain.copy(deep=True)
        self.output.attrs = attrs['_global']
        # Number of timesteps
        if self.disagg:
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
        # Input preprocessing
        in_preprocess = {"ascii" : self.vic_in_preprocess,
                         "binary" : self.vic_in_preprocess,
                         "netcdf" : self.netcdf_in_preprocess}
        in_preprocess[MetSim.params['in_format']](job_list)

        # Output preprocessing
        out_preprocess = {"ascii": self.ascii_out_preprocess,
                          "netcdf": self.netcdf_out_preprocess}
        out_preprocess[MetSim.params['out_format']]()
        self.ready = True

    def launch(self, job_list):
        """Farm out the jobs to separate processes"""
        nprocs = MetSim.params['nprocs']
        self.pool = Pool(processes=nprocs)

        # Split the input into chunks to run in parallel
        locations = np.array_split(list(zip(self.i_idx, self.j_idx)), 
                nprocs * cnst.CHUNK_SIZE)

        # Do the forcing generation and disaggregation if required
        status = []
        for loc_chunk in locations:
            stat = self.pool.apply_async(
                wrap_run,
                args=(self.method.run, loc_chunk, self.met_data, self.disagg),
                callback=self._unpack_results)
            status.append(stat)
        self.pool.close()

        # Check that everything worked
        [stat.get() for stat in status]
        self.pool.join()

    def _unpack_results(self, results):
        """Put results into the master dataset"""
        for result in results:
            locs, df = result
            i = locs['lat']
            j = locs['lon']
            for varname in self.params['out_vars']:
                self.output[varname].values[:, i, j] = df[varname].values

    def run(self, locations):
        """
        Kicks off the disaggregation and queues up data for IO
        """
        results = []
        for i, j in locations:
            print("Processing {} {}".format(i, j))
            locs = dict(lat=i, lon=j)
            ds = self.met_data.isel(lat=i, lon=j)  # fix this indexing
            lat = ds['lat'].values
            elev = ds['elev'].values
            df = ds.drop(['lat', 'lon', 'elev']).to_dataframe()
            df = self.method.run(df, MetSim.params, elev=elev,
                    lat=lat, disagg=self.disagg)
            results.append((locs, df))
        self._unpack_results(results)

    def find_elevation(self, lat: float, lon: float) -> float:
        """ Use the domain file to get the elevation """
        i = int(np.abs(self.lat - lat).argmin())
        j = int(np.abs(self.lon - lon).argmin())
        return self.elev.isel(lat=i, lon=j)

    def _validate_setup(self):
        """Updates the global parameters dictionary"""
        errs = [""]

        # Make sure there's some input
        if not len(self.params['forcing']):
            errs.append("Requires input forcings to be specified")

        # Parameters that can't be empty strings or None
        non_empty = ['method', 'domain', 'out_dir',
                     'start', 'stop', 'time_step', 'out_format',
                     'in_format', 't_max_lr', 't_min_lr']
        for each in non_empty:
            if self.params[each] is None or self.params[each] == '':
                errs.append("Cannot have empty value for {}".format(each))

        # Check for required input variable specification
        required_in = ['t_min', 't_max', 'prec']
        for each in required_in:
            if each not in self.met_data.variables:
                errs.append("Input requires {}".format(each))

        # Convert data types as necessary
        self.params['t_max_lr'] = float(self.params['t_max_lr'])
        self.params['t_min_lr'] = float(self.params['t_min_lr'])

        # Make sure that we are going to write out some data
        if not len(self.params['out_vars']):
            errs.append("Output variable list must not be empty")

        # Check that the parameters specified are available
        opts = {'mtclim_swe_corr': [True, False],
                'lw_cloud': ['default', 'cloud_deardorff'],
                'lw_type': ['default', 'tva', 'anderson',
                            'brutsaert', 'satterlund',
                            'idso', 'prata']}
        for k, v in opts.items():
            if not MetSim.params[k] in v:
                errs.append("Invalid option given for {}".format(k))

        # If any errors, raise and give a summary
        if len(errs) > 1:
            raise Exception("\n  ".join(errs))

    def netcdf_in_preprocess(self, filename):
        """Get the extent and load the data"""
        self.met_data = self.read(filename)

    def vic_in_preprocess(self, job_list):
        """Process all files to find spatial extent"""
        # Creates the master dataset which will be used to parallelize
        self.met_data = xr.Dataset(coords={'time' : MetSim.params['dates'],
                                      'lon' : self.lon,
                                      'lat' : self.lat},
                              attrs={'n_days' : len(MetSim.params['dates'])})
        shape = (len(MetSim.params['dates']), len(self.lat), len(self.lon))

        self.met_data['elev'] = (('lat', 'lon'),
                           np.full((len(self.lat), len(self.lon)), np.nan))
        for var in MetSim.params['in_vars']:
            self.met_data[var] = (('time', 'lat', 'lon'),np.full(shape, np.nan))

        # Fill in the data
        for job in job_list:
            _, lat, lon = os.path.basename(job).split("_")
            i = np.unique(self.i_idx)[list(self.lat.values).index(float(lat))]
            j = np.unique(self.j_idx)[list(self.lon.values).index(float(lon))]
            ds = self.read(job)
            self.met_data['elev'].values[i, j] = self.elev.values[i, j]
            for var in MetSim.params['in_vars']:
                self.met_data[var].values[:, i, j] = ds[var].values

    def ascii_out_preprocess(self):
        """Dummy function"""
        pass

    def netcdf_out_preprocess(self):
        """Initialize the output file"""
        print("Initializing netcdf...")
        self.output_filename = os.path.join(self.params['out_dir'], "forcing.nc")

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

    def write_ascii(self):
        """Write out as ASCII to the output file"""
        print("Writing ascii...")
        shape = self.output.dims
        for i, j in itertools.product(range(shape['lat']), range(shape['lon'])):
            if self.output.mask[i, j] > 0:
                lat = self.output.lat.values[i]
                lon = self.output.lon.values[j]
                fname = os.path.join(self.params['out_dir'], 
                            "forcing_{}_{}.csv".format(lat, lon))
                self.output.isel(lat=i, lon=j)[self.params[
                    'out_vars']].to_dataframe().to_csv(fname)
        
    def read(self, fpath: str) -> xr.Dataset:
        """
        Dispatch to the right function based on the file extension
        """
        dispatch = {
                'binary': self.read_binary,
                'netcdf': self.read_netcdf,
                'ascii': self.read_ascii
                }
        return dispatch[MetSim.params['in_format']](fpath)

    def read_binary(self, fpath: str) -> xr.Dataset:
        """ Reads a binary forcing file (VIC 4 format) """
        dates = pd.date_range(MetSim.params['start'], MetSim.params['stop'])
        n_days = len(dates)
        type_lookup = {'signed': 'h', 'unsigned': 'H'}
        # Pack these for nicer syntax in the loop
        var_names = MetSim.params['in_vars'].keys()
        data_list = [[] for var in MetSim.params['in_vars'].keys()]
        n_vars = len(var_names)
        params = [p for p in MetSim.params['in_vars'].values()]
        scales = [float(s.split()[0]) for s in params]
        datatypes = [type_lookup[s.split()[-1]] for s in params]
        with open(fpath, 'rb') as f:
            i = 0
            points_read = 0
            points_needed = 4*n_days
            while points_read != points_needed:
                bytes = f.read(2)
                if bytes:
                    # Get correct variable and data type with i,
                    # then unpack & scale
                    data_list[i].append(
                        struct.unpack(datatypes[i], bytes)[0] / scales[i])
                    i = (i + 1) % n_vars
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

        # Assemble the dataset
        data_dict = {c[0]: (['time'], c[1]) for c in zip(var_names, data_list)}
        data_dict['elev'] = (['lon', 'lat'], params['elev'])
        data_dict['day_of_year'] = (['time'], dates.dayofyear)
        df = xr.Dataset(data_dict,
                        coords={'lon': [params['lon']],
                                'lat': [params['lat']],
                                'time': dates},
                        attrs={'n_days': params['n_days']})
        return df

    def read_ascii(self, fpath:str) -> xr.Dataset:
        """Read in an ascii forcing file"""
        dates = pd.date_range(MetSim.params['start'], MetSim.params['stop'])
        ds = pd.read_table(fpath, header=None, delim_whitespace=True,
                names=MetSim.params['in_vars'].keys()).head(len(dates))
        ds.index = dates
        return ds

    def read_netcdf(self, fpath: str) -> xr.Dataset:
        """
        Read in a NetCDF file and add elevation information
        """
        dates = pd.date_range(MetSim.params['start'], MetSim.params['stop'])
        n_days = len(dates)
        ds = xr.open_dataset(fpath)
        ds = ds.sel(time=slice(MetSim.params['start'], MetSim.params['stop']))
        ds.rename(MetSim.params['in_vars'], inplace=True)

        # Add elevation and day of year data
        ds['elev'] = xr.Variable(['lon', 'lat'], self.elev.values)
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
        ds = met_data.isel(lat=i, lon=j)  # fix this indexing
        lat = ds['lat'].values
        elev = ds['elev'].values
        df = ds.drop(['lat', 'lon', 'elev']).to_dataframe()
        df = func(df, MetSim.params, elev=elev, lat=lat, disagg=disagg)
        results.append((locs, df))
    return results

