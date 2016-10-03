"""
Options for MetSim
"""

import os
import sys
import argparse
import metsim

def parse(args):
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--config', 
            default=None)

    parser.add_argument('-n', '--n-processes',
            default=1)

    return parser.parse_args()


def init(opts):
    """
    TODO
    """
    print(opts)
    if not os.path.isfile(opts.config):
        exit("Invalid configuration given.  Use `ms -h` for more information.")
    metsim.config = opts.config
    metsim.proc_count = opts.n_processes    

