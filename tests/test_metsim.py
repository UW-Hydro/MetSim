#!/usr/bin/env python3
"""
Unit tests for MetSim
"""

import os
import pytest
import numpy as np
import pandas as pd
import xarray as xr

from metsim.metsim import MetSim
import metsim.constants as const

# Parameters to test over
in_fmts = ['ascii', 'binary', 'netcdf']
out_fmts = ['ascii', 'netcdf']
methods = ['mtclim']

# Where datasets for each input type are found
data_locations = {'netcdf' : './tests/data/test.nc',
                  'ascii' : './tests/data/ascii/',
                  'binary' : './tests/data/binary/'}

# Dates to run over
dates = {'netcdf' : (pd.datetime(1950,1,1), pd.datetime(1950,1,31)),
         'binary' : (pd.datetime(1949,1,1), pd.datetime(1949,12,31)),
         'ascii' : (pd.datetime(1949,1,1), pd.datetime(2005,12,31))}

# All values should be in these ranges
data_ranges = {'temp' : (-50,40),
               'prec' : (0,8),
               'shortwave' : (0,1000),
               'longwave' : (0,450),
               'wind' : (0,10),
               'vapor_pressure' : (0,1.4),
               'rel_humid' : (0,100)}



@pytest.fixture(params=in_fmts)
def in_format(request):
    """Input formats - see in_fmts`"""
    return request.param


@pytest.fixture(params=out_fmts)
def out_format(request):
    """Output formats - see `out_fmts`"""
    return request.param


@pytest.fixture(params=methods)
def method(request):
    """Generation methods - see `methods`"""
    return request.param


@pytest.fixture()
def domain_file():
    """Domain file containing elevation data"""
    return "./tests/data/domain.nc"


@pytest.fixture()
def test_params(in_format, out_format, method):
    """Assemble the parameters for each combo"""
    start = dates[in_format][0]
    stop = dates[in_format][1]
    in_vars = ['prec', 't_max', 't_min', 'wind']
    out_dir = "./tmp"
    lr = 0.0065
    params = {'start' : start,
              'stop' : stop,
              'in_vars' : in_vars,
              'in_format' : in_format,
              'out_format' : out_format,
              'domain' : './tests/data/domain.nc',
              'forcings' : data_locations[in_format],
              'method' : method,
              'time_step' : "60",
              't_max_lr' : lr,
              't_min_lr' : lr,
              'out_dir' : out_dir
              }
    return params


@pytest.fixture()
def test_setup(test_params, domain_file):
    """Tests the setup of the MetSim object"""
    in_fmt = test_params['in_format']
    loc = data_locations[in_fmt]

    # Get the files and make sure the right amount exist
    if in_fmt == 'binary':
        data_files = [os.path.join(loc, f) for f in os.listdir(loc)]
        assert len(data_files) == 16
    elif in_fmt == 'ascii':
        data_files = [os.path.join(loc, f) for f in os.listdir(loc)]
        assert len(data_files) == 1
    else:
        data_files = loc
        assert data_files == './tests/data/test.nc'

    # Test construction - should not yet be ready to run
    ms = MetSim(test_params)
    assert not ms.ready

    # Load data and ensure the ready flag has been set
    ms.load(data_files, domain_file)
    assert ms.ready

    # Check to see that the data is valid
    assert type(ms.met_data) is xr.Dataset

    return ms


def test_mtclim(test_setup):
    """Tests the ability to run successfully"""
    # Here we only test a single grid cell
    loc = test_setup.locations[0]
    strloc = "{}_{}".format(loc[0], loc[1])
    n_days = len(test_setup.met_data.time)

    # Run the forcing generation, but not the disaggregation
    daily = test_setup.run([loc], disagg=False)[strloc]
    assert len(daily) == n_days
    for var in ['prec', 't_max', 't_min', 't_day', 'wind', 'elev',
                'dayl', 'swrad', 'tskc', 'pet', 'vapor_pressure']:
        assert var in daily

    # Now test the disaggregation as well as forcing generation
    hourly = test_setup.run([loc])[strloc]
    assert len(hourly) == (n_days * const.HOURS_PER_DAY)+1
    for var in test_setup.params['out_vars']:
        assert var in hourly
        l, h = data_ranges[var]
        assert hourly[var].between(l, h).all()


def test_disaggregation_values():
    """Tests to make sure values are being generated correctly"""
    # Set parameters 
    params = {'start' : dates['binary'][0],
              'stop' : dates['binary'][1],
              'in_vars' : ['prec', 't_max', 't_min', 'wind'],
              'in_format' : 'binary',
              'out_format' : 'ascii',
              'domain' : './tests/data/domain.nc',
              'forcings' : ['./tests/data/binary/data_48.3125_-120.5625'],
              'method' : 'mtclim',
              'time_step' : "60",
              't_max_lr' : 0.00065,
              't_min_lr' : 0.00065,
              'out_dir' : "./tmp" 
              }   
    # The location we will test against
    loc = (48.3125,-120.5625)

    # Set up the MetSim object
    ms = MetSim(params)
    ms.load(params['forcings'], params['domain'])

    # Run MetSim and load in the validated data
    out = ms.run([loc])["{}_{}".format(loc[0], loc[1])]
    good = pd.read_table('./tests/data/validated_48.3125_-120.5625', index_col=0)

    # Make sure the data comes out right
    assert type(out) is pd.DataFrame
    for var in out.keys():
        # Check to make sure each variable is within  
        # 0.005 validated data
        assert np.allclose(out[var], good[var], 1e-5, 0.005)

