"""
Handles the synchronization of multiple processes for MetSim
"""

import os
import struct
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset

from metsim import configuration

class MetSim(object):
    """
    MetSim handles the distribution of jobs that write to a common file
    by launching muliple processes and queueing up their writeback so that 
    work can be done while IO is happening.
    """

    def __init__(self, domain:str, params:dict):
        """
        Constructor
        """
        # Record parameters and broadcast to the config module
        self.params = params
        configuration.update(params) 

        # Get the necessary information from the domain 
        domain = Dataset(configuration.PARAMS['domain'], 'r')
        self.elev = np.array(domain['elev'])
        self.lats = np.array(domain['lat'])
        self.lons = np.array(domain['lon'])
        domain.close()

        # If we are outputting netcdf, intiialize the file
        if params.get('out_format', 'netcdf') == 'netcdf':
            self.init_netcdf()


    def run(self, method, forcings):
        """
        Kicks off the disaggregation and queues up data for IO
        """
        # Where we will store the results for writing out
        out_dict = {}
           
        # Coerce the forcings into a list if it's not one
        if type(forcings) is not list:
            forcings = [forcings]
           
        # Do the forcing generation and dissaggregation if required
        for forcing in forcings:
            met_data = self.read(forcing)
            for i in range(len(met_data.lat)):
                out_dict[i] = {}
                for j in range(len(met_data.lon)):
                    out_dict[i][j] = method.run(
                                met_data.isel(lat=[i], lon=[j])
                                .to_dataframe(), disagg=True)
                    # TODO: Remove this
                    print(i,j) 
            # Write out
            self.write(out_dict)

        if self.params.get('out_format', 'netcdf') == 'netcdf':
            self.output.close()
           

    def find_elevation(self, lat: float, lon: float) -> float:
        """ Use the domain file to get the elevation """
        lat_idx = np.abs(self.lats - lat).argmin()
        lon_idx = np.abs(self.lons - lon).argmin()            
        return self.elev[lat_idx, lon_idx]

   
    def write(self, data: dict):
        """
        Dispatch to the right function based on the configuration given 
        """
        dispatch = {
                'netcdf' : self.write_netcdf,
                'ascii'  : self.write_ascii
                }
        dispatch[self.params.get('out_format', 'netcdf').lower()](data)
         

    def init_netcdf(self):
        """Initialize the output file"""
        print("Initializing netcdf...")
        out_file = os.path.join(self.params['out_dir'], "forcing.nc")
        self.output = Dataset(out_file, 'w')

        # Number of timesteps 
        n_ts = len(pd.date_range(self.params['start'], 
                self.params['stop'] + pd.Timedelta('1 days'), 
                freq=self.params['time_step']+'T'))
        
        # Create dimensions
        self.output.createDimension('lat', len(self.lats))
        self.output.createDimension('lon', len(self.lons))
        self.output.createDimension('time', n_ts) 

        # Create output variables
        for varname in self.params['out_vars']:
            self.output.createVariable(varname, 'd', ('time', 'lat','lon'))


    def write_netcdf(self, data:dict):
        """Write out as NetCDF to the output file"""
        print("Writing netcdf...")
        lats = data.keys()
        for lat in lats:
            lons = data[lat].keys()
            for lon in lons:
                print(lat, lon)
                for varname in self.params['out_vars']:
                    self.output.variables[varname][:, lat, lon] = data[lat][lon][varname].values


    def write_ascii(self, data: dict):
        """Write out as ASCII to the output file"""
        print("Writing ascii...")
        out_base = os.path.join(self.params['out_dir'], "forcing")
        lats = data.keys()
        for lat in lats:
            lons = data[lat].keys()
            for lon in lons:
                print(lat, lon)
                data[lat][lon].to_csv('_'.join([out_base, str(lat), str(lon)]), sep='\t')
        

    def read_binary(self, fpath: str) -> xr.Dataset:
        """ Reads a binary forcing file (VIC 4 format) """
        dates = pd.date_range(self.params['start'], self.params['stop'])
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
        configuration.update(params)
        params['elev'] = [[self.find_elevation(params['lat'], params['lon'])]]
        df = xr.Dataset({"precip"      : (['time'], precip), 
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


    def read_netcdf(self, fpath:str) -> xr.Dataset:
        """
        Read in a NetCDF file and add elevation information 
        """
        # TODO: FIXME: Finish this
        dates = pd.date_range(self.params['start'], self.params['stop'])
        n_days = len(dates)
        ds = xr.open_dataset(fpath)
        ds = ds.sel(time=slice(self.params['start'], self.params['stop']))
        ds.rename({"Prec" : "precip",
                   "Tmax" : "t_max",
                   "Tmin" : "t_min"}, inplace=True)

        # Add elevation and day of year data
        ds['elev'] = (['lat','lon'], self.elev)
        ds['day_of_year'] = (['time'], dates.dayofyear)

        # Update the configuration
        configuration.update({"n_days" : n_days})

        #TODO: If the domain file and input file don't have the same grid
        #      structure this will have to be implemented
        #in_lats, in_lons = np.meshgrid(self.lats, self.lons)
        #out_lats, out_lons = np.meshgrid(ds.lat, ds.lon)
        #self.elev = basemap.interp(self.elev, self.lats, self.lons, out_lats, out_lons, order=1)
        return ds 
    
    
    def read(self, fpath:str) -> xr.Dataset:
        """
        Dispatch to the right function based on the file extension 
        """
        ext_to_fun = {
                '.bin'   : self.read_binary,
                '.nc'    : self.read_netcdf,
                '.nc4'   : self.read_netcdf
                }
        return ext_to_fun.get(os.path.splitext(fpath)[-1], self.read_binary)(fpath) 
    

