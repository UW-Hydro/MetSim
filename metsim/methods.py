"""
Contains functions for disaggregation of meteorological data
"""

import os
import time

import metsim
import metsim.io

def mtclim(flist, writable=True):
    for f in flist:
        data = f
        metsim.io.sync_io(writable, metsim.io.hold_lock, data, "na")


# A mapping from the config variable to the function handle
mapping = {
           'mtclim' : mtclim
          }

