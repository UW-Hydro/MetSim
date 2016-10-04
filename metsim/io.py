"""
Handles IO for MetSim
"""

import os
import struct
from configparser import ConfigParser

def read_config(fname):
    """
    TODO
    """
    if not os.path.isfile(fname):
        return dict()
    cfp = ConfigParser()
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
    vars = [precip, t_min, t_max, wind]

    # Data types referred to: 'H' - unsigned short ; 'h' - short
    types = ['H', 'h', 'h', 'h']
    with open(fname, 'rb') as f:
        i = 0
        while True:
            bytes = f.read(2)
            if bytes:
                # Get correct variable and data type with i, then unpack
                vars[i].append(struct.unpack(types[i], bytes)[0])
                i = (i+1)%4
            else:
                break
    return vars


def write_to_netcdf(data, fname):
    """
    TODO
    """
    pass


def init_netcdf(fname):
    """
    TODO
    """
    pass

