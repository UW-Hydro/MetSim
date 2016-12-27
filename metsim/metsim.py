"""
Handles the synchronization of multiple processes for MetSim
"""

import os
import struct
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
from mpl_toolkits import basemap 
from multiprocessing import Value, Process

from metsim import configuration

class MetSim(object):
    """
    MetSim handles the distribution of jobs that write to a common file
    by launching muliple processes and queueing up their writeback so that 
    work can be done while IO is happening.
    """

    def __init__(self, domain, params):
        """
        Constructor
        """
        # Builds the infrastructure to keep track of jobs
        self.writable = Value('b', True, lock=False)

        self.params = params
        configuration.update(params) 

        # Keep a handle to the domain file
        domain = Dataset(configuration.PARAMS['domain'], 'r')
        self.elev = np.array(domain['elev'])
        self.lats = np.array(domain['lat'])
        self.lons = np.array(domain['lon'])
        domain.close()

        # Set up the distribution of jobs and create process handles
        #self.jobs = [job_list[i:i+job_size] for i in range(0, n_jobs, job_size)]
        #self.process_handles = [
        #         Process(target=self.run, args=(self.method, job_list))
        #         for job_list in self.jobs
        #        ]
        self.process_handles = []


    def run(self, method, forcings):
        """
        Kicks off the disaggregation and queues up data for IO
        """
        if type(forcings) is not list:
            forcings = [forcings]
        for forcing in forcings:
            met_data = self.read(forcing)
            for i in range(len(met_data.lat)):
                for j in range(len(met_data.lon)):
                    # FIXME: This needs to be fixed.
                    #        See little notebook for details
                    temp = method.run(
                                met_data.isel(lat=[i], lon=[j])
                                .to_dataframe(), disagg=True)
                    
                    if j == 100: exit()
                    print(i, j)

            
    def launch_processes(self):
        """ Launches all processes built in the constructor """
        for p in self.process_handles:
            p.start()
        for p in self.process_handles:
            p.join()


    def find_elevation(self, lat: float, lon: float) -> float:
        """ Use the domain file to get the elevation """
        lat_idx = np.abs(self.lats - lat).argmin()
        lon_idx = np.abs(self.lons - lon).argmin()            
        return self.elev[lat_idx, lon_idx]

   
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


    def read_netcdf(self, fpath, n_days=-1) -> xr.Dataset:
        """
        TODO
        """
        # TODO: FIXME: Finish this
        dates = pd.date_range(self.params['start'], self.params['stop'])
        n_days = len(dates)
        ds = xr.open_dataset(fpath)
        ds = ds.sel(time=slice(self.params['start'], self.params['stop']))
        ds.rename({"Prec" : "precip",
                   "Tmax" : "t_max",
                   "Tmin" : "t_min"}, inplace=True)
        ds['elev'] = (['lat','lon'], self.elev)
        ds['day_of_year'] = (['time'], dates.dayofyear)
        configuration.update({"n_days" : n_days})
        #in_lats, in_lons = np.meshgrid(self.lats, self.lons)
        #out_lats, out_lons = np.meshgrid(ds.lat, ds.lon)
        #self.elev = basemap.interp(self.elev, self.lats, self.lons, out_lats, out_lons, order=1)
        return ds 
    
    
    def read(self, fpath) -> xr.Dataset:
        """
        Dispatch to the right function based on the file extension 
        """
        ext_to_fun = {
                '.bin'   : self.read_binary,
                '.nc'    : self.read_netcdf,
                '.nc4'   : self.read_netcdf
                }
        return ext_to_fun.get(os.path.splitext(fpath)[-1], self.read_binary)(fpath) 
    

