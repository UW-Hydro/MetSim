"""
Options for MetSim
"""

import os
import sys
import glob
import argparse
import metsim
from metsim import io
from metsim import methods

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
    metsim.config = io.read_config(opts.config)

    #NOTE: This will silently override invalid methods in the configuration file
    metsim.method = methods.mapping.get(metsim.config['Output']['disagg_method'], 'mtclim')

    # Split up the list of forcing files for the desired number of processes
    # Note here that we auto force the number of processes to be at most
    # the number of forcing files.
    metsim.forcing_files = os.listdir(metsim.config['IO']['forcing_dir']) 
    metsim.proc_count = min(opts.n_processes, len(metsim.forcing_files))    
    chunk = int(len(metsim.forcing_files) / metsim.proc_count)
    metsim.forcing_chunks = [metsim.forcing_files[i:i+chunk] for i in 
            range(0,len(metsim.forcing_files),chunk)]

