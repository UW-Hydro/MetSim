"""
Handles IO for MetSim
"""

import os
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
    pass


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

