#!/usr/bin/env python3
"""
Unit tests for MetSim
"""

import os
import pytest
import pandas as pd
import xarray as xr

from metsim.metsim import MetSim
import metsim.constants as const

in_fmts = ['ascii', 'binary', 'netcdf']
out_fmts = ['ascii', 'netcdf']

data_locations = {'netcdf' : './tests/data/test.nc',
                  'ascii'  : './tests/data/ascii/',
                  'binary' : './tests/data/binary/'}

# TODO: Fix these for each dataset
dates = {'netcdf' : (pd.datetime(1950,1,1), pd.datetime(1950,1,31)),
         'binary' : (pd.datetime(1949,1,1), pd.datetime(1949,12,31)),
         'ascii' : (pd.datetime(1949,1,1), pd.datetime(2005,12,31))}

methods = ['mtclim']


@pytest.fixture(params=in_fmts)
def in_format(request):
    return request.param


@pytest.fixture(params=out_fmts)
def out_format(request):
    return request.param


@pytest.fixture(params=methods)
def method(request):
    return request.param


@pytest.fixture()
def domain_file():
    return "./tests/data/domain.nc"


@pytest.fixture()
def test_params(in_format, out_format, method):
    start = dates[in_format][0] 
    stop = dates[in_format][1]
    in_vars = ['prec', 'wind', 't_max', 't_min']
    out_dir = "./tmp"
    lr = 0.00065
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


