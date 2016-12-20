"""
Handles the synchronization of multiple processes for MetSim
"""

import os
import pandas as pd
from multiprocessing import Value, Process

from metsim import io
from metsim import configuration

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
            forcing = io.read(job, len(dates))
            forcing.set_dates(dates)
            forcing.generate_met_forcings(method) 
            forcing.disaggregate()
            io.sync_io(io.write_ascii, forcing.met_data, self.writable, 
                       os.path.join(os.path.basename(job))) 


    def launch_processes(self):
        """ Launches all processes built in the constructor """
        for p in self.process_handles:
            p.start()
        for p in self.process_handles:
            p.join()

