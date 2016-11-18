"""
Handles IO for MetSim
"""

import os
import csv
import time
import struct
import pandas as pd
from configparser import ConfigParser
from metsim.forcing import Forcing

def read_config(fpath: str, n_days=-1) -> ConfigParser:
    """
    Reads the configuration file given.  If the
    file does not exist, an empty ConfigParser
    will be returned.
    """
    cfp = ConfigParser()
    if os.path.isfile(fpath):
        cfp.read(fpath)
    return cfp 
    

def read_ascii(fpath, n_days=-1) -> Forcing:
    """
    TODO
    """
    precip = []
    t_min = []
    t_max = []
    wind = []

    with open(fpath, 'r') as f:
        for line in csv.reader(f, delimiter='\t'):
            print(line)
            break
    # TODO: FIXME: Finish this
    return Forcing(None, None)


def read_binary(fpath: str, n_days=-1) -> Forcing:
    """ Reads a binary forcing file (VIC 4 format) """
    precip = [] # Short unsigned int
    t_max  = [] # Short int
    t_min  = [] # Short int
    wind   = [] # Short int
    
    # Pack these for nicer syntax in the loop
    var_name = [precip, t_max, t_min, wind]
    scale = [40.0, 100.0, 100.0, 100.0]

    # Data types referred to: 'H' - unsigned short ; 'h' - short
    types = ['H', 'h', 'h', 'h']
    with open(fpath, 'rb') as f:
        i = 0
        points_read = 0
        points_needed = 4*n_days
        while points_read != points_needed:
            bytes = f.read(2)
            if bytes:
                # Get correct variable and data type with i, then unpack & scale
                var_name[i].append(struct.unpack(types[i], bytes)[0]/scale[i])
                i = (i+1)%4
                points_read += 1
            else:
                break

    # Binary forcing files have naming format $NAME_$LAT_$LON
    param_list = os.path.basename(fpath).split(".")[-1].split("_")
    params = {"name"   : param_list[0], 
              "lat"    : param_list[1], 
              "lon"    : param_list[2],
              "n_days" : n_days}
    df = pd.DataFrame(data={"precip" : precip, 
                            "t_min"  : t_min, 
                            "t_max"  : t_max, 
                            "wind"   : wind})
    return Forcing(df, params) 



def read_netcdf(fpath, n_days=-1) -> Forcing:
    """
    TODO
    """
    # TODO: FIXME: Finish this
    return Forcing(None, None) 


def read(fpath, n_days=-1) -> Forcing:
    """
    Dispatch to the right function based on the file extension 
    """
    ext_to_fun = {
            '.txt'   : read_ascii,
            '.ascii' : read_ascii,
            '.bin'   : read_binary,
            ''       : read_binary,
            '.nc'    : read_netcdf,
            '.nc4'   : read_netcdf,
            '.conf'  : read_config,
            '.ini'   : read_config
            }
    return ext_to_fun[os.path.splitext(fpath)[-1]](fpath, n_days) 


def write(data, fpath):
    return fpath.split('.')[-1]


def write_ascii(data, fpath):
    """
    TODO
    """
    data.to_csv(fpath, sep='\t')


def write_netcdf(data, fpath):
    """
    TODO
    """
    hold_lock(data, fpath) 


def init_netcdf(fpath):
    """
    TODO
    """
    pass


def hold_lock(data, fpath, timeout=10):
    """
    A dummy method to hold the IO lock for some amount of time.
    """ 
    print("Holding lock")
    time.sleep(timeout)
    print(data)


def sync_io(io_function, io_data, writable, io_fpath):
    """
    Simple wrapper to make it easy to check when to do IO
    """
    # Wait until "lock" is released
    if writable.value == 0:
        print("Waiting for lock...")
        time.sleep(3)
    # Grab the lock and start doing IO
    writable.value = False
    io_function(io_data, io_fpath)
    writable.value = True


