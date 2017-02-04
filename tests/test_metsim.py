#!/usr/bin/env python3
"""
Unit tests for MetSim
"""

import os
import pytest
import pandas as pd

from metsim.metsim import MetSim

in_fmts = ['ascii', 'binary', 'netcdf']
out_fmts = ['ascii', 'netcdf']

data_locations = {'netcdf' : './tests/data/test.nc',
                  'ascii'  : './tests/data/ascii/',
                  'binary' : './tests/data/binary/'}

# TODO: Fix these for each dataset
dates = {'netcdf' : (pd.datetime(2000,1,1), pd.datetime(2001,1,1)),
         'binary' : (pd.datetime(2000,1,1), pd.datetime(2001,1,1)),
         'ascii' : (pd.datetime(2000,1,1), pd.datetime(2001,1,1))}

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
              't_max_lr' : lr,
              't_min_lr' : lr,
              'out_dir' : out_dir
              }
    return params


def test_setup(test_params, domain_file):
    in_fmt = test_params['in_format']
    loc = data_locations[in_fmt]
    if in_format == 'ascii' or in_format == 'binary':
        data_files = os.listdir(loc) 
        assert len(data_files) == 16
    else:
        data_files = loc
        assert data_files == './tests/data/test.nc'

    ms = MetSim(test_params)
    assert not ms.ready
    ms.load(data_files, domain_file)
    assert ms.ready


