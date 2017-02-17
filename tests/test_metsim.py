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

# Domain files to use
domain_files = {'netcdf' : './tests/data/domain.nc',
                'ascii' : './tests/data/stehekin.nc',
                'binary' : './tests/data/stehekin.nc'}

# Dates to run over
dates = {'netcdf' : (pd.datetime(1950,1,1), pd.datetime(1950,1,31)),
         'binary' : (pd.datetime(1949,1,1), pd.datetime(1949,12,31)),
         'ascii' : (pd.datetime(1949,1,1), pd.datetime(1949,12,31))}

# Domain vars
domain_section = {'netcdf': {'lat': 'lat', 
                             'lon': 'lon', 
                             'mask': 'mask', 
                             'elev': 'elev'},
                  'binary': {'lat': 'lat', 
                             'lon': 'lon', 
                             'mask': 'mask', 
                             'elev': 'elev'},
                  'ascii': {'lat': 'lat', 
                            'lon': 'lon', 
                            'mask': 'mask', 
                            'elev': 'elev'}
                  }

# Input vars
in_vars_section = {'netcdf': {'Prec': 'prec', 
                              'Tmax': 't_max', 
                              'Tmin': 't_min'},
                   'binary': {'prec': '40.0 unsigned',
                              't_max': '100.0 signed',
                              't_min': '100.0 signed',
                              'wind': '100.0 signed'},
                   'ascii': {'prec': 'prec', 
                             't_max': 't_max', 
                             't_min': 't_min', 
                             'wind': 'wind'}}


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
    in_vars = in_vars_section[in_format] 
    domain_vars = domain_section[in_format]
    out_dir = "./tmp"
    lr = 0.0065
    params = {'start' : start,
              'stop' : stop,
              'in_vars' : in_vars,
              'in_format' : in_format,
              'out_format' : out_format,
              'domain' : domain_files[in_format],
              'method' : method,
              'time_step' : "60",
              't_max_lr' : lr,
              't_min_lr' : lr,
              'out_dir' : out_dir,
              'in_vars' : in_vars,
              'domain_vars' : domain_vars
              }
    return params


@pytest.fixture()
def test_setup(test_params, domain_file):
    """Tests the setup of the MetSim object"""
    in_fmt = test_params['in_format']
    loc = data_locations[in_fmt]

    # Get the files and make sure the right amount exist
    if in_fmt == 'binary' or in_fmt == 'ascii':
        data_files = [os.path.join(loc, f) for f in os.listdir(loc)]
        assert len(data_files) == 16
    else:
        data_files = loc
        assert data_files == './tests/data/test.nc'
    test_params['forcings'] = data_files

    # Test construction - should not yet be ready to run
    ms = MetSim(test_params)
    assert not ms.ready

    return ms


def test_mtclim(test_setup):
    """Tests the ability to run successfully"""
    # Here we only test a single grid cell
    data_files = test_setup.params['forcings']
    daily_out_vars = ['prec', 't_max', 't_min', 't_day', 
                      'wind', 'dayl', 'swrad', 
                      'tskc', 'pet', 'vapor_pressure']
    hourly_out_vars = ['prec', 'temp', 'shortwave', 'longwave',
                      'vapor_pressure', 'rel_humid']
    
    # Load data and ensure the ready flag has been set
    test_setup.params['time_step'] = 1440
    test_setup.params['out_vars'] = daily_out_vars
    test_setup.load(data_files)
    loc = test_setup.locations[0]
    assert test_setup.ready

    # Check to see that the data is valid
    assert type(test_setup.met_data) is xr.Dataset
    n_days = len(test_setup.met_data.time)

    # Run the forcing generation, but not the disaggregation
    test_setup.run([loc])
    daily = test_setup.output
    assert type(test_setup.output) is xr.Dataset
    assert len(daily.time) == n_days
    for var in daily_out_vars:
        assert var in daily
    
     
    # Now test the disaggregation as well as forcing generation
    test_setup.params['time_step'] = 60
    test_setup.params['out_vars'] = hourly_out_vars 
    test_setup.load(data_files)
    assert test_setup.ready

    # Check to see that the data is valid
    assert type(test_setup.met_data) is xr.Dataset
    
    test_setup.run([loc])
    hourly = test_setup.output.isel(lat=loc[0], lon=loc[1]).to_dataframe()
    assert len(hourly) == (n_days * const.HOURS_PER_DAY)+1
    for var in test_setup.params['out_vars']:
        assert var in hourly
        l, h = data_ranges[var]
        assert hourly[var].between(l, h).all()


def test_disaggregation_values():
    """Tests to make sure values are being generated correctly"""
    # Set parameters 
    loc = data_locations['binary']
    data_files = [os.path.join(loc, f) for f in os.listdir(loc)]
    out_vars = ['prec', 'temp', 'shortwave', 'longwave',
                'vapor_pressure', 'wind', 'rel_humid']
    params = {'start' : dates['binary'][0],
              'stop' : dates['binary'][1],
              'in_format' : 'binary',
              'out_format' : 'ascii',
              'domain' : './tests/data/stehekin.nc',
              'forcings' : data_files, 
              'method' : 'mtclim',
              'time_step' : "60",
              't_max_lr' : 0.00065,
              't_min_lr' : 0.00065,
              'out_dir' : "./tmp",
              'out_vars' : out_vars,
              'in_vars' : in_vars_section['binary'],
              'domain_vars': domain_section['binary']
              }   
    # The location we will test against
    loc = (1,4)

    # Set up the MetSim object
    ms = MetSim(params)
    ms.load(params['forcings'])

    # Run MetSim and load in the validated data
    ms.run([loc])
    out = ms.output.isel(lat=loc[0], lon=loc[1]).to_dataframe()[:-1][out_vars]
    good = pd.read_table('./tests/data/validated_48.3125_-120.5625',
            names=out_vars)
    good.index = out.index

    # Make sure the data comes out right
    assert type(out) is pd.DataFrame
    for var in ms.params['out_vars']:
        # Check to make sure each variable has normalized  
        # rmse of less than 0.02 
        h = max([good[var].max(), out[var].max()])
        l = min([good[var].min(), out[var].min()])
        nrmse = np.sqrt((good[var] - out[var]).pow(2).mean())/(h-l)
        assert nrmse < 0.02
