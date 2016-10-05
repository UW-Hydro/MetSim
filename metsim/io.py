"""
Handles IO for MetSim
"""

import os
import time
import struct
from configparser import ConfigParser

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
    pass


def read_binary_forcng(fname):
    """
    TODO
    """
    #TODO: There are scale factors that still need to be applied
    precip = [] # Short unsigned int
    t_min  = [] # Short int
    t_max  = [] # Short int
    wind   = [] # Short int
    
    # Pack these for nicer syntax in the loop
    var_name = [precip, t_min, t_max, wind]
    scale = [40.0, 100.0, 100.0, 100.0]

    # Data types referred to: 'H' - unsigned short ; 'h' - short
    types = ['H', 'h', 'h', 'h']
    with open(fname, 'rb') as f:
        i = 0
        while True:
            bytes = f.read(2)
            if bytes:
                # Get correct variable and data type with i, then unpack
                var_name[i].append(struct.unpack(types[i], bytes)[0]/scale[i])
                i = (i+1)%4
            else:
                break
    return vars


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


