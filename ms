#!/usr/bin/env python3

"""
Command line tool for MetSim
"""

import os
import sys
import json
import argparse
import pandas as pd
from configparser import ConfigParser


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
    config = dict(config['MetSim'])
    out_dir = config.get('out_dir', '')
    method = config.get('method', '')
    
    # If the forcing variable is a directory, scan it for files
    if os.path.isdir(config.get('forcing', '')):
        forcing_files = [os.path.join(config.get('forcing', ''),fn) for fn in
                         next(os.walk(config.get('forcing', '')))[2]]
    else:
        forcing_files = config.get('forcing', '')

    # We assume there is only one domain file
    domain_file = config.get('domain', None)

    # Generate the date range that will be put into the data frame
    start_date, start_hour = config.get('start', '').split(':')
    start_year, start_month, start_day = start_date.split('/')
    end_year, end_month, end_day = config.get('stop', '').split('/')
    start = pd.datetime(int(start_year),
                        int(start_month),
                        int(start_day),
                        int(start_hour))
    stop  = pd.datetime(int(end_year),
                        int(end_month),
                        int(end_day))
    to_list = lambda s : json.loads(s.replace("'", '"'))
    
    config.update({"stop" : stop,
                   "start" : start,
                   "nprocs" : opts.n_processes,
                   "method" : method,
                   "out_dir" : out_dir,
                   "domain" : domain_file,
                   "forcings" : forcing_files})
    config['in_vars'] = to_list(config.get('in_vars', '[]')) 
    config['out_vars'] = to_list(config.get('out_vars', '[]')) 
    config['t_min_lr'] = config.get('t_min_lr', "")
    config['t_max_lr'] = config.get('t_max_lr', "")
    config = {k : v for k, v in config.items() if v != [] and v is not None}
    return config 


def main():
    """Runs MetSim"""
    from metsim.metsim import MetSim
    setup = init(parse(sys.argv[1:]))
    ms = MetSim(setup)
    ms.launch(setup['forcings'])


if __name__ == '__main__':
    main()

