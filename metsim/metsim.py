"""
Handles the synchronization of multiple processes for MetSim
"""

import time
from multiprocessing import Value, Process

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
        self.method = analysis_method
        n_jobs = len(job_list) 
        job_size = int(n_jobs / min(n_processes, n_jobs))
        self.jobs = [job_list[i:i+job_size] for i in range(0, n_jobs, job_size)]
        self.process_handles = [
                 Process(target=self.method, args=(job_list, self.writable))
                 for job_list in self.jobs
                ]

    def run(self):
        for p in self.process_handles:
            p.start()




