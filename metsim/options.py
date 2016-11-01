"""
Options for MetSim
"""

import os
import sys
import glob
import argparse
import pandas as pd
import metsim
from metsim import io
from metsim import defaults

def parse(args):
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--config', 
            default=None)

    parser.add_argument('-n', '--n-processes',
            default=1, type=int)

    return parser.parse_args()


def init(opts):
    """
    TODO
    """
    if not os.path.isfile(opts.config):
        exit("Invalid configuration given.  Use `ms -h` for more information.")
    metsim.config = io.read(opts.config)
    metsim.input_format = metsim.config['IO']['force_format']
    metsim.out_dir = metsim.config['IO']['out_dir']
    try:
        os.mkdir(metsim.out_dir)
    except:
        print("ERROR: Could not create directory: " + metsim.out_dir)
        exit()

    #NOTE: This will silently override invalid methods in the configuration file
    metsim.method = defaults.METHODS.get(metsim.config['Output']['disagg_method'], 
                                        'mtclim')

    # Split up the list of forcing files for the desired number of processes
    # Note here that we auto force the number of processes to be at most
    # the number of forcing files.
    metsim.forcing_files = [os.path.join(metsim.config['IO']['forcing_dir'], f) for f in 
        os.listdir(metsim.config['IO']['forcing_dir'])] 
    metsim.proc_count = min(opts.n_processes, len(metsim.forcing_files))    
    chunk = int(len(metsim.forcing_files) / metsim.proc_count)
    metsim.forcing_chunks = [metsim.forcing_files[i:i+chunk] for i in 
            range(0,len(metsim.forcing_files),chunk)]

    # Generate the date range that will be put into the data frame
    metsim.start = pd.datetime(int(metsim.config['IO']['start_year']),
                        int(metsim.config['IO']['start_month']),
                        int(metsim.config['IO']['start_day']),
                        int(metsim.config['IO']['start_hour']))
    metsim.stop  = pd.datetime(int(metsim.config['IO']['end_year']),
                        int(metsim.config['IO']['end_month']),
                        int(metsim.config['IO']['end_day']))

