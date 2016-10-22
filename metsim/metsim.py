"""
Handles the synchronization of multiple processes for MetSim
"""

import os
import time
import numpy as np
import pandas as pd
from multiprocessing import Value, Process

import metsim.io

class MetSim(object):
    """
    MetSim handles the distribution of jobs that write to a common file
    by launching muliple processes and queueing up their writeback so that 
    work can be done while IO is happening.
    """

    def __init__(self, analysis_method, job_list, n_processes):
        """
        Constructor
        """
        # Builds the infrastructure to keep track of jobs
        self.writable = Value('b', True, lock=False)
        
        # Set up the distribution of jobs and create process handles
        self.method = analysis_method.run
        n_jobs = len(job_list) 
        job_size = int(n_jobs / min(n_processes, n_jobs))
        self.jobs = [job_list[i:i+job_size] for i in range(0, n_jobs, job_size)]
        self.process_handles = [
                 Process(target=self.run, args=(self.method, job_list))
                 for job_list in self.jobs
                ]


    def run(self, method, job_list):
        """
        Kicks off the disaggregation and queues up data for IO
        """
        data = pd.DataFrame()
        for job in job_list:
            #TODO: Fix this so the correct file format is chosen - multiple dispatch?
            forcing = metsim.io.read(job)
            data = pd.concat([data, self.method(forcing)])
        #FIXME: This needs to either return data or set a class variable
        self.disaggregate(data)
        metsim.io.sync_io(metsim.io.write_ascii, forcing, self.writable, 
                    os.path.join(metsim.out_dir, os.path.basename(job))) 


    def launch_processes(self):
        """ Launches all processes built in the constructor """
        for p in self.process_handles:
            p.start()
        for p in self.process_handles:
            p.join()


    def disaggregate(self, df_daily):
        """
        TODO
        """
        #FIXME: Can we use multiple dispatch here as well?
        dates_hourly = pd.date_range(metsim.start, metsim.stop, freq='H') 
        df_hourly = pd.DataFrame(index=dates_hourly)
        df_hourly = (pd.concat(self._disagg_temp(   df_daily, df_hourly))
                       .concat(self._disagg_precip( df_daily, df_hourly))
                       .concat(self._disagg_thermal(df_daily, df_hourly))
                       .concat(self._disagg_wind(   df_daily, df_hourly)))
    
 
    def _disagg_precip(self, df_daily, df_hourly):
        pass


    def _disagg_temp(self, df_daily, df_hourly):
        pass


    def _disagg_thermal(self, df_daily, df_hourly):
        pass


    def _disagg_wind(self, df_daily, df_hourly):
        pass


   

