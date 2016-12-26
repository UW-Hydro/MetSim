"""
Handles the synchronization of multiple processes for MetSim
"""

import os
import struct
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from multiprocessing import Value, Process

from metsim import configuration
from metsim.forcing import Forcing

class MetSim(object):
    """
    MetSim handles the distribution of jobs that write to a common file
    by launching muliple processes and queueing up their writeback so that 
    work can be done while IO is happening.
    """

    def __init__(self, forcings, domain, method, params):
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
        self.run(method, forcings)
        #self.jobs = [job_list[i:i+job_size] for i in range(0, n_jobs, job_size)]
        #self.process_handles = [
        #         Process(target=self.run, args=(self.method, job_list))
        #         for job_list in self.jobs
        #        ]
        self.process_handles = []


    def run(self, method, job_list):
        """
        Kicks off the disaggregation and queues up data for IO
        """
        for job in job_list:
            dates = pd.date_range(self.params['start'], self.params['stop'])
            forcing = self.read(job, len(dates))
            forcing.set_dates(dates)
            forcing.generate_met_forcings(method) 
            forcing.disaggregate()

            
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

   
    def read_binary(self, fpath: str, n_days=-1) -> Forcing:
        """ Reads a binary forcing file (VIC 4 format) """
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
        params['elev'] = self.find_elevation(params['lat'], params['lon'])
        df = pd.DataFrame(data={"precip" : precip, 
                                "t_min"  : t_min, 
                                "t_max"  : t_max, 
                                "wind"   : wind})
        return Forcing(df, params) 


    def read_netcdf(self, fpath, n_days=-1) -> Forcing:
        """
        TODO
        """
        # TODO: FIXME: Finish this
    
        return Forcing(None, None) 
    
    
    def read(self, fpath, n_days=-1) -> Forcing:
        """
        Dispatch to the right function based on the file extension 
        """
        ext_to_fun = {
                '.bin'   : self.read_binary,
                '.nc'    : self.read_netcdf,
                '.nc4'   : self.read_netcdf
                }
        return ext_to_fun.get(os.path.splitext(fpath)[-1], self.read_binary)(fpath, n_days) 
    

