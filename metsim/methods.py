"""
Contains functions for disaggregation of meteorological data
"""

import os
import time

import metsim
import metsim.io

def mtclim(fname):
    # Read input forcing
    # Do the disagg
    # Return a pandas or xray object
    pass

# A mapping from the config variable to the function handle
mapping = {
           'mtclim' : mtclim
          }

