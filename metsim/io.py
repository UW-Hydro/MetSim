"""
Handles IO for MetSim
"""

import os
import csv
import time
import struct
import pandas as pd
from configparser import ConfigParser
import metsim

def read_config(fname):
    """
    TODO
    """
    cfp = ConfigParser()
    if os.path.isfile(fname):
        cfp.read(fname)
    return cfp 
    

def read_ascii_forcing(fname):
    """
    TODO
    """
    precip = []
    t_min = []
    t_max = []
    wind = []

    with open(fname, 'r') as f:
        for line in csv.reader(f, delimiter='\t'):
            print(line)
            break


def read_binary_forcng(fname):
    """
    TODO
    """
    precip = [] # Short unsigned int
    t_max  = [] # Short int
    t_min  = [] # Short int
    wind   = [] # Short int
    
    # Pack these for nicer syntax in the loop
    var_name = [precip, t_max, t_min, wind]
    scale = [40.0, 100.0, 100.0, 100.0]

    # Data types referred to: 'H' - unsigned short ; 'h' - short
    types = ['H', 'h', 'h', 'h']
    with open(fname, 'rb') as f:
        i = 0
        points_read = 0
        points_needed = 4*len(metsim.dates)
        while points_read != points_needed:
            bytes = f.read(2)
            if bytes:
                # Get correct variable and data type with i, then unpack
                var_name[i].append(struct.unpack(types[i], bytes)[0]/scale[i])
                i = (i+1)%4
                points_read += 1
            else:
                break
    df = pd.DataFrame(data={"precip" : precip, 
                            "t_min"  : t_min, 
                            "t_max"  : t_max, 
                            "wind"   : wind},
                      index=metsim.dates)
    return df

def read_netcdf_forcing(fname):
    """
    TODO
    """
    pass


def write_ascii(data, fname):
    """
    TODO
    """
    data.to_csv(fname, sep='\t')


def write_netcdf(data, fname):
    """
    TODO
    """
    hold_lock(data, fname) 


def init_netcdf(fname):
    """
    TODO
    """
    pass


def hold_lock(data, fname, timeout=10):
    """
    A dummy method to hold the IO lock for some amount of time.
    """ 
    print("Holding lock")
    time.sleep(timeout)
    print(data)


def sync_io(io_function, io_data, writable, io_fname):
    """
    Simple wrapper to make it easy to check when to do IO
    """
    # Wait until "lock" is released
    if writable.value == 0:
        print("Waiting for lock...")
        time.sleep(3)
    # Grab the lock and start doing IO
    writable.value = False
    io_function(io_data, io_fname)
    writable.value = True


