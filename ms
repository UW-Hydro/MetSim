#!/usr/bin/env python3

"""
Command line tool for MetSim
"""

import os
import sys
import argparse
import pandas as pd
from configparser import ConfigParser

from metsim.configuration import METHODS as methods

def parse(args):
    """Parse the command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', default=None)
    parser.add_argument('-n', '--n-processes', default=1, type=int)
    return parser.parse_args()


def init(opts):
    """Initialize some information based on the options & config"""
    if not os.path.isfile(opts.config):
        exit("Invalid configuration given.  Use `ms -h` for more information.")
    config = ConfigParser()
    config.read(opts.config)
    config = config['MetSim']
    out_dir = config['out_dir']
    
    try:
        os.mkdir(out_dir)
    except:
        print("ERROR: Could not create directory: " + out_dir)
        exit(1)

    #NOTE: This will silently override invalid methods in the configuration file
    method = methods.get(config['method'], 'mtclim')
   
    # If the forcing variable is a directory, scan it for files
    if os.path.isdir(config['forcing']):
        forcing_files = [os.path.join(config['forcing'],fn) for fn in 
                         next(os.walk(config['forcing']))[2]]
    else:
        forcing_files = [config['forcing']]
    
    # We assume there is only one domain file
    domain_file = config['domain']

    # Generate the date range that will be put into the data frame
    start_date, start_hour = config['start'].split(':')
    start_year, start_month, start_day = start_date.split('/')
    end_year, end_month, end_day = config['stop'].split('/')
    start = pd.datetime(int(start_year),
                        int(start_month),
                        int(start_day),
                        int(start_hour))
    stop  = pd.datetime(int(end_year),
                        int(end_month),
                        int(end_day))
    setup = {"stop"      : stop,
             "start"     : start, 
             "method"    : method, 
             "domain"    : domain_file,
             "forcings"  : forcing_files,
             "in_vars"   : config['in_vars'],
             "out_vars"  : config['out_vars'],
             "time_step" : config['time_step']}
    return setup 


def main():
    """Runs MetSim"""
    from metsim.metsim import MetSim 
    setup = init(parse(sys.argv[1:])) 
    ms = MetSim(setup['forcings'], setup['domain'], setup['method'], setup)
    ms.launch_processes()


if __name__ == '__main__':
    main()
    
