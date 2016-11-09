"""
Handles the synchronization of multiple processes for MetSim
"""

import os
import time
import numpy as np
import pandas as pd
from multiprocessing import Value, Process

import metsim.io
from metsim.disaggregate import disaggregate

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
        self.run(self.method, job_list)
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
            dates = pd.date_range(metsim.start, metsim.stop)
            forcing = metsim.io.read(job, len(dates))
            forcing = forcing.set_index(dates)
            forcing['day_of_year'] = dates.dayofyear
            metsim.n_days = len(forcing['day_of_year'])
            forcing = self.method(forcing)
        # Discard the daily data in favor of hourly data, then write
        forcing = disaggregate(forcing)
        metsim.io.sync_io(metsim.io.write_ascii, forcing, self.writable, 
                    os.path.join(metsim.out_dir, os.path.basename(job))) 


    def launch_processes(self):
        """ Launches all processes built in the constructor """
        for p in self.process_handles:
            p.start()
        for p in self.process_handles:
            p.join()

