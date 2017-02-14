"""
Handles the synchronization of multiple processes for MetSim
"""

import os
import struct
import itertools
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
from multiprocessing import Process, Manager

from metsim.methods import mtclim

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
        "method":'',
        "domain":'',
        "out_dir":'',
        "start":'',
        "stop":'',
        "time_step":'',
        "out_format":'',
        "in_format" : None,
        "base_elev" : 0,
        "t_max_lr": '',
        "t_min_lr": '',
        "out_vars" : ['temp', 'prec', 'shortwave', 'longwave',
                      'wind', 'vapor_pressure', 'rel_humid']
    }

    missing_data_msg = """No data to run on!  Use the `load`
                          method to set the data to operate on.
                          For more information about how to use
                          MetSim refer to the documentation."""

    def __init__(self, params:dict):
        """
        Constructor
        """
        # Record parameters
        self.update(params)
        MetSim.params['dates'] = pd.date_range(params['start'], params['stop'])
        self.output = None
        self.met_data = None
        self.ready = False

        # Create the data structures for writeout
        self.manager = Manager()
        self.queue = self.manager.Queue()


    def load(self, job_list, domain):
        """Preprocess the input data"""
        # Get the necessary information from the domain
        self.domain = Dataset(domain, 'r')
        self.domain_lats = np.array(self.domain['lat'])
        self.domain_lons = np.array(self.domain['lon'])
        self.elev = np.array(self.domain['elev'])

        # Used to dispatch to the correct preprocessing functions
        in_preprocess = {"ascii" : self.ascii_in_preprocess,
                         "binary" : self.binary_in_preprocess,
                         "netcdf" : self.netdf_in_preprocess}
        self.met_data = in_preprocess[MetSim.params['in_format']](job_list)

        self.ready = True


    def launch(self, job_list):
        """Farm out the jobs to separate processes"""
        # Load in the data
        self.load(job_list, self.params['domain'])

        # Start up the IO Process
        io_process = Process(target=self._create_io_process, args=(self.queue,))
        io_process.start()

        # Split the jobs up into groups based on the number
        # of desired processes
        self.locations = np.array_split(np.array(self.locations),
                                        MetSim.params['nprocs'])
        process_handles = [Process(target=self.run, args=(locs, self.queue))
                           for locs in self.locations]

        # Runs everything
        for p in process_handles:
            p.start()
        for p in process_handles:
            p.join()

        # Close everything out
        self.queue.put("done")
        io_process.join()


    def run(self, locations, queue=None, disagg=True):
        """
        Kicks off the disaggregation and queues up data for IO
        """
        # Where we will store the results for writing out
        out_dict = {}
        method = MetSim.methods[MetSim.params.get('method', 'mtclim')]

        # Do the forcing generation and disaggregation if required
        for i, j in locations:
            print("Processing {} {}".format(i, j))
            out_dict["{}_{}".format(i,j)]  = (
                    method.run(self.met_data.sel(lat=[i], lon=[j]).to_dataframe(),
                               MetSim.params, disagg=disagg))

        # Output to the queue if we have an IO process
        if queue is not None:
            queue.put(out_dict)

        return out_dict


    def find_elevation(self, lat: float, lon: float) -> float:
        """ Use the domain file to get the elevation """
        lat_idx = np.abs(self.domain_lats - lat).argmin()
        lon_idx = np.abs(self.domain_lons - lon).argmin()
        return self.domain['elev'][lat_idx, lon_idx]


    def update(self, new_params):
        """Updates the global parameters dictionary"""
        MetSim.params.update(new_params)


    def netdf_in_preprocess(self, filename):
        """Get the extent and load the data"""
        # Get the information for splitting up the job
        in_forcing = Dataset(filename, 'r')
        self.locations = list(itertools.product(in_forcing['lat'], in_forcing['lon']))
        self.lats = np.array(in_forcing['lat'])
        self.lons = np.array(in_forcing['lon'])
        self.lat_idx = {str(lat): i for i, lat in enumerate(self.lats)}
        self.lon_idx = {str(lon): j for j, lon in enumerate(self.lons)}
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


    def ascii_out_preprocess(self):
        """Dummy function"""
        pass


    def binary_out_preprocess(self):
        """Dummy function"""
        pass


    def netcdf_out_preprocess(self):
        """Initialize the output file"""
        print("Initializing netcdf...")
        out_dir = MetSim.params['out_dir']
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        out_file = os.path.join(out_dir, "forcing.nc")
        self.output = Dataset(out_file, 'w')

        # Number of timesteps
        n_ts = len(pd.date_range(MetSim.params['start'],
                   MetSim.params['stop'] + pd.Timedelta('1 days'),
                   freq="{}T".format(MetSim.params['time_step'])))

        # Create dimensions
        self.output.createDimension('lat', len(self.lats))
        self.output.createDimension('lon', len(self.lons))
        self.output.createDimension('time', n_ts)

        # Create output variables
        for varname in MetSim.params['out_vars']:
            self.output.createVariable(varname, 'd', ('time', 'lat','lon'))


    def _create_io_process(self, queue):
        """Launch the IO process"""
        out_preprocess = {"ascii" : self.ascii_out_preprocess,
                          "binary" : self.binary_out_preprocess,
                          "netcdf" : self.netcdf_out_preprocess}
        out_preprocess[MetSim.params['out_format']]()
        # Wait for data to come in
        while 1:
            data = queue.get()
            if data == "done":
                break
            self.write(data)

        # Close the dataset when done
        if self.output:
            self.output.close()


    def _validate_parameters(self):
        """
        Make sure that all of the necessary parameters
        have been provided.
        """
        validation_msg = ""
        raise Exception(validation_msg)


    def write(self, data: dict):
        """
        Dispatch to the right function based on the configuration given
        """
        dispatch = {
                'netcdf' : self.write_netcdf,
                'ascii'  : self.write_ascii
                }
        dispatch[MetSim.params.get('out_format', 'netcdf').lower()](data)


    def write_netcdf(self, data:dict):
        """Write out as NetCDF to the output file"""
        print("Writing netcdf...")
        for l in data.keys():
            lat, lon = l.split("_")
            for varname in MetSim.params['out_vars']:
                i, j = self.lat_idx[lat], self.lon_idx[lon]
                self.output.variables[varname][:, i, j] = (
                        data[l][varname].values)


    def write_ascii(self, data: dict):
        """Write out as ASCII to the output file"""
        print("Writing ascii...")
        out_dir = MetSim.params.get('out_dir', './results/')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        for l in data.keys():
            lat, lon = l.split("_")
            out_file = '_'.join(["forcing", str(lat), str(lon)])
            data[l].to_csv(os.path.join(out_dir, out_file), sep='\t')


    def read(self, fpath:str) -> xr.Dataset:
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
        precip = [] # Short unsigned int
        t_max  = [] # Short int
        t_min  = [] # Short int
        wind   = [] # Short int

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
                    # Get correct variable and data type with i, then unpack & scale
                    var_name[i].append(struct.unpack(types[i], bytes)[0]/scale[i])
                    i = (i+1)%4
                    points_read += 1
                else:
                    break

        # Binary forcing files have naming format $NAME_$LAT_$LON
        param_list = os.path.basename(fpath).split("_")
        params = {"name"   : param_list[0],
                  "lat"    : float(param_list[1]),
                  "lon"    : float(param_list[2]),
                  "n_days" : int(n_days)}
        MetSim.params.update(params)
        params['elev'] = [[self.find_elevation(params['lat'], params['lon'])]]
        df = xr.Dataset({"prec"      : (['time'], precip),
                         "t_min"       : (['time'], t_min),
                         "t_max"       : (['time'], t_max),
                         "wind"        : (['time'], wind),
                         "elev"        : (['lon', 'lat'], params['elev']),
                         "day_of_year" : (['time'], dates.dayofyear)},
                        coords={'lon'  : [params['lon']],
                                'lat'  : [params['lat']],
                                'time' : dates},
                        attrs={'n_days' : params['n_days']})
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
        ds.rename({"Prec" : "prec",
                   "Tmax" : "t_max",
                   "Tmin" : "t_min"}, inplace=True)

        # Add elevation and day of year data
        ds['elev'] = (['lat','lon'], self.elev)
        ds['day_of_year'] = (['time'], dates.dayofyear)

        # Update the configuration
        MetSim.params.update({"n_days" : n_days})
        return ds


