#!/usr/bin/env python3
"""
Unit tests for MetSim
"""

import os
import subprocess
import tempfile
from collections import OrderedDict

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from datetime import datetime

import metsim.cli.ms as cli
import metsim.metsim
from metsim.metsim import MetSim
from metsim import io


class DummyOpts:
    def __init__(self, config):
        self.config = config
        self.scheduler = 'threading'
        self.verbose = False
        self.num_workers = 1


# Parameters to test over
in_fmts = ['ascii', 'binary', 'netcdf']
methods = ['mtclim']
timesteps = [1440, 30]

# Where datasets for each input type are found
data_locations = {'netcdf': './metsim/data/test.nc',
                  'ascii': './metsim/data/ascii/',
                  'binary': './metsim/data/binary/'}

# Domain files to use
domain_files = {'netcdf': './metsim/data/domain.nc',
                'ascii': './metsim/data/stehekin.nc',
                'binary': './metsim/data/stehekin.nc'}

# State files to use
state_files = {'netcdf': './metsim/data/state_nc.nc',
               'ascii': './metsim/data/state_vic.nc',
               'binary': './metsim/data/state_vic.nc'}

# Dates to run over
dates = {'netcdf': (datetime(1950, 1, 1), datetime(1950, 1, 31)),
         'binary': (datetime(1949, 1, 1), datetime(1949, 12, 31)),
         'ascii': (datetime(1949, 1, 1), datetime(1949, 12, 31))}

# Domain vars
domain_section = {'netcdf': OrderedDict(lat='lat', lon='lon', mask='mask',
                                        elev='elev'),
                  'binary': OrderedDict(lat='lat', lon='lon', mask='mask',
                                        elev='elev'),
                  'ascii': OrderedDict(lat='lat', lon='lon', mask='mask',
                                       elev='elev')}

# Input vars
in_vars_section = {'netcdf': OrderedDict(Prec='prec', Tmax='t_max',
                                         Tmin='t_min', wind='wind'),
                   'binary': OrderedDict([('prec', '40.0 unsigned'),
                                          ('t_max', '100.0 signed'),
                                          ('t_min', '100.0 signed'),
                                          ('wind', '100.0 signed')]),
                   'ascii': OrderedDict([('prec', 'prec'),
                                         ('t_max', 't_max'),
                                         ('t_min', 't_min'),
                                         ('wind', 'wind')])}

# All values should be in these ranges
data_ranges = {'temp': (-50, 40),
               'prec': (0, 8),
               'shortwave': (0, 1000),
               'longwave': (0, 450),
               'wind': (0, 10),
               'vapor_pressure': (0, 2),
               'air_pressure': (0, 101325),
               'spec_humid': (0, 2),
               'rel_humid': (0, 100)}


@pytest.fixture(params=in_fmts)
def in_format(request):
    """Input formats - see in_fmts`"""
    return request.param


@pytest.fixture(params=methods)
def method(request):
    """Generation methods - see `methods`"""
    return request.param


@pytest.fixture(params=timesteps)
def timestep(request):
    """Generation methods - see `methods`"""
    return request.param



@pytest.fixture()
def domain_file():
    """Domain file containing elevation data"""
    return "./tests/data/domain.nc"


@pytest.fixture()
def test_params(in_format, method, timestep):
    """Assemble the parameters for each combo"""
    start = dates[in_format][0]
    stop = dates[in_format][1]
    in_vars = in_vars_section[in_format]
    domain_vars = domain_section[in_format]
    out_dir = tempfile.mkdtemp('results')
    out_prefix = "forcing"
    params = {'start': start,
              'stop': stop,
              'in_vars': in_vars,
              'forcing_fmt': in_format,
              'domain_fmt': 'netcdf',
              'state_fmt': 'netcdf',
              'domain': domain_files[in_format],
              'state': state_files[in_format],
              'method': method,
              'calender': 'standard',
              'time_step': timestep,
              'time_grouper': None,
              'out_dir': out_dir,
              'out_state': os.path.join(out_dir, 'state.nc'),
              'out_prefix': out_prefix,
              'forcing_vars': in_vars,
              'domain_vars': domain_vars}
    return params


@pytest.fixture()
def test_setup(test_params, domain_file):
    """Tests the setup of the MetSim object"""
    in_fmt = test_params['forcing_fmt']
    loc = data_locations[in_fmt]

    # Get the files and make sure the right amount exist
    if in_fmt == 'binary' or in_fmt == 'ascii':
        data_files = [os.path.join(loc, f) for f in os.listdir(loc)]
        assert len(data_files) == 16
    else:
        data_files = loc
        assert data_files == './metsim/data/test.nc'
    test_params['forcing'] = data_files

    # Test construction
    ms = MetSim(test_params)
    return ms


#def test_mtclim(test_setup):
#    """Tests the ability to run successfully"""
#    # Here we only test a single grid cell
#    daily_out_vars = ['prec', 't_max', 't_min', 'wind', 'shortwave',
#                      'tskc', 'pet', 'vapor_pressure']
#    hourly_out_vars = ['prec', 'temp', 'shortwave', 'longwave',
#                       'vapor_pressure', 'wind', 'rel_humid', 'spec_humid',
#                       'air_pressure']
#
#    # Load data and ensure the ready flag has been set
#    test_setup.params['time_step'] = 1440
#    test_setup.params['out_vars'] = daily_out_vars
#
#    # Check to see that the data is valid
#    assert isinstance(test_setup.met_data, xr.Dataset)
#    n_days = len(test_setup.met_data.time)
#
#    # Run the forcing generation, but not the disaggregation
#    test_setup.run()
#    daily = test_setup.open_output()
#    print(daily)
#    assert isinstance(daily, xr.Dataset)
#    assert len(daily.time) == n_days
#    for var in daily_out_vars:
#        assert var in daily
#
#    # Now test the disaggregation as well as forcing generation
#    test_setup.params['time_step'] = 60
#    test_setup.params['out_vars'] = hourly_out_vars
#
#    # Check to see that the data is valid
#    assert isinstance(test_setup.met_data, xr.Dataset)
#
#    test_setup.run()
#    hourly = test_setup.open_output().isel(lat=2, lon=2).to_dataframe()
#    assert len(hourly) == (n_days * const.HOURS_PER_DAY)
#    for var in test_setup.params['out_vars']:
#        assert var in hourly
#        l, h = data_ranges[var]
#        vl = min(hourly[var].values)
#        vh = max(hourly[var].values)
#        print(var, vl, vh, l, h)
#        assert hourly[var].between(l, h).all()
#
#    # Now test sub-hourly disaggregation
#    test_setup.params['time_step'] = 30
#    test_setup.run()
#    half_hourly = test_setup.open_output().isel(lat=1, lon=3).to_dataframe()
#    assert len(half_hourly) == (2 * n_days * const.HOURS_PER_DAY)

def test_time_offset():
    """Tests to make sure that the time_offset option works"""
    loc = data_locations['binary']
    data_files = [os.path.join(loc, f) for f in os.listdir(loc)]
    out_vars = ['prec', 'temp', 'shortwave', 'longwave', 'vapor_pressure',
                'wind', 'rel_humid', 'spec_humid', 'air_pressure']
    out_dir = '.'
    params = {'start': dates['binary'][0],
              'stop': dates['binary'][1],
              'forcing_fmt': 'binary',
              'domain_fmt': 'netcdf',
              'state_fmt': 'netcdf',
              'domain': './metsim/data/stehekin.nc',
              'state': './metsim/data/state_vic.nc',
              'forcing': data_files,
              'method': 'mtclim',
              'scheduler': 'threading',
              'time_step': "60",
              'out_dir': out_dir,
              'out_state': os.path.join(out_dir, 'state.nc'),
              'out_vars': {n: metsim.metsim.available_outputs[n]
                           for n in out_vars},
              'forcing_vars': in_vars_section['binary'],
              'domain_vars': domain_section['binary']
              }
    params1 = dict()
    params1.update(params)
    params2 = dict()
    params2.update(params)
    params1['period_ending'] = False
    params2['period_ending'] = True

    # Set up the MetSim object
    ms1 = MetSim(params1)
    ms2 = MetSim(params2)
    assert ms1._times[1:] == ms2._times[:-1]


def test_variable_rename():
    """Tests to make sure that variable renaming works"""
    loc = data_locations['binary']
    data_files = [os.path.join(loc, f) for f in os.listdir(loc)]
    out_dir = '.'
    params = {'start': dates['binary'][0],
              'stop': dates['binary'][1],
              'forcing_fmt': 'binary',
              'domain_fmt': 'netcdf',
              'state_fmt': 'netcdf',
              'domain': './metsim/data/stehekin.nc',
              'state': './metsim/data/state_vic.nc',
              'forcing': data_files,
              'method': 'mtclim',
              'scheduler': 'threading',
              'time_step': "60",
              'out_dir': out_dir,
              'out_state': os.path.join(out_dir, 'state.nc'),
              'out_vars': {
                  'prec': {'out_name': 'pptrate'},
                  'shortwave': {'out_name': 'SWRadAtm'}},
              'forcing_vars': in_vars_section['binary'],
              'domain_vars': domain_section['binary']
              }
    ms = MetSim(params)
    ms.run()
    ds = ms.open_output()
    assert 'pptrate' in ds.variables
    assert 'SWRadAtm' in ds.variables


def test_passthrough():
    """Test to make sure passing through previously estimated
       variables doesn't alter the values from disaggregation"""
    loc = data_locations['binary']
    data_files = [os.path.join(loc, f) for f in os.listdir(loc)]
    out_dir = '.'
    params = {'start': dates['binary'][0],
              'stop': dates['binary'][1],
              'forcing_fmt': 'binary',
              'domain_fmt': 'netcdf',
              'state_fmt': 'netcdf',
              'domain': './metsim/data/stehekin.nc',
              'state': './metsim/data/state_vic.nc',
              'forcing': data_files,
              'method': 'mtclim',
              'scheduler': 'threading',
              'time_step': "60",
              'out_dir': out_dir,
              'out_state': os.path.join(out_dir, 'state.nc'),
              'out_vars': {'prec': {'out_name': 'prec'},
                           'temp': {'out_name': 'temp'},
                           'longwave': {'out_name': 'longwave'},
                           'vapor_pressure': {'out_name': 'vapor_pressure'},
                           'shortwave': {'out_name': 'shortwave'}},
              'forcing_vars': in_vars_section['binary'],
              'domain_vars': domain_section['binary']}

    # Second run will be to use mtclim and hourly disagg
    params1 = dict()
    params1.update(params)
    params1['out_prefix'] = 'mtclim'
    ms1 = MetSim(params1)
    ms1.run()
    with ms1.open_output() as ds:
        mtclim_ds = ds.load()

    # Third run will be to use passthrough and hourly disagg
    # with input data from teh first run
    params2 = dict()
    params2.update(params)
    params2['method'] = 'passthrough'
    params2['out_prefix'] = 'passthrough'
    params2['forcing_vars'] = OrderedDict(
        prec='prec', t_max='t_max', t_min='t_min',
        wind='wind', shortwave='shortwave', vapor_pressure='vapor_pressure')
    params2['forcing_fmt'] = 'netcdf'
    params2['forcing'] = './metsim/data/passthrough.nc'
    ms2 = MetSim(params2)
    ms2.run()
    with ms2.open_output() as ds:
        passthrough_ds = ds.load()

    tol = 1e-4
    assert np.allclose(passthrough_ds['shortwave'].mean(),
                       mtclim_ds['shortwave'].mean(), atol=tol)
    assert np.allclose(passthrough_ds['vapor_pressure'].mean(),
                       mtclim_ds['vapor_pressure'].mean(), atol=tol)
    assert np.allclose(passthrough_ds['longwave'].mean(),
                       mtclim_ds['longwave'].mean(), atol=tol)



def test_unit_conversion():
    """Tests to make sure that variable renaming works"""
    loc = data_locations['binary']
    data_files = [os.path.join(loc, f) for f in os.listdir(loc)]
    out_dir = '.'
    params = {'start': dates['binary'][0],
              'stop': dates['binary'][1],
              'forcing_fmt': 'binary',
              'domain_fmt': 'netcdf',
              'state_fmt': 'netcdf',
              'domain': './metsim/data/stehekin.nc',
              'state': './metsim/data/state_vic.nc',
              'forcing': data_files,
              'method': 'mtclim',
              'scheduler': 'threading',
              'time_step': "60",
              'out_dir': out_dir,
              'out_state': os.path.join(out_dir, 'state.nc'),
              'out_vars': {
                  'prec': {'out_name': 'pptrate',
                           'units': 'mm s-1'},
                  'temp': {'out_name': 'airtemp',
                           'units': 'K'}},
              'forcing_vars': in_vars_section['binary'],
              'domain_vars': domain_section['binary']}

    params1 = dict()
    params1.update(params)
    params2 = dict()
    params2.update(params)
    params2['out_vars'] = {
        'prec': {'out_name': 'pptrate',
                 'units': 'mm timestep-1'},
        'temp': {'out_name': 'airtemp',
                 'units': 'C'}}
    ms1 = MetSim(params1)
    ms1.run()
    ds1 = ms1.open_output().load()
    ds1.close()
    time_step = int(params['time_step'])
    sec_per_min = 60.
    tol = 1e-4

    ms2 = MetSim(params2)
    ms2.run()
    ds2 = ms2.open_output().load()

    assert np.allclose(ds1['airtemp'].mean(),
                       ds2['airtemp'].mean()+273.15, atol=tol)
    assert np.allclose(time_step * sec_per_min * ds1['pptrate'].mean(),
                       ds2['pptrate'].mean(), atol=tol)



def test_disaggregation_values():
    """Tests to make sure values are being generated correctly"""
    # Set parameters
    loc = data_locations['binary']
    data_files = [os.path.join(loc, f) for f in os.listdir(loc)]
    out_vars = ['prec', 'temp', 'shortwave', 'longwave', 'vapor_pressure',
                'wind', 'rel_humid', 'spec_humid', 'air_pressure']
    out_dir = tempfile.mkdtemp('results')
    params = {'start': dates['binary'][0],
              'stop': dates['binary'][1],
              'forcing_fmt': 'binary',
              'domain_fmt': 'netcdf',
              'state_fmt': 'netcdf',
              'domain': './metsim/data/stehekin.nc',
              'state': './metsim/data/state_vic.nc',
              'forcing': data_files,
              'method': 'mtclim',
              'scheduler': 'threading',
              'time_step': "60",
              'out_dir': out_dir,
              'out_state': os.path.join(out_dir, 'state.nc'),
              'out_vars': {n: metsim.metsim.available_outputs[n]
                           for n in out_vars},
              'forcing_vars': in_vars_section['binary'],
              'domain_vars': domain_section['binary']
              }
    # The location we will test against
    loc = (1, 4)

    def check_data(out, good, tol=0.03):
        assert isinstance(out, pd.DataFrame)
        for var in ms.params['out_vars'].keys():
            # Check to make sure each variable has normalized
            # rmse of less than 0.02
            h = max([good[var].max(), out[var].max()])
            l = min([good[var].min(), out[var].min()])
            nrmse = np.sqrt((good[var] - out[var]).pow(2).mean()) / (h - l)
            print(var, nrmse)
            assert nrmse < tol

    # Set up the MetSim object
    ms = MetSim(params)

    # Run MetSim and load in the validated data
    ms.run()
    ds = ms.open_output()
    out = ds.isel(lat=loc[0], lon=loc[1]).to_dataframe()[out_vars]
    good = pd.read_table('./metsim/data/validated_48.3125_-120.5625',
                         names=out_vars)
    good.index = out.index

    # Make sure the data comes out right
    check_data(out, good)
    ds.close()

    # Now do 3 hourly
    params['time_step'] = '180'
    ms = MetSim(params)
    ms.run()
    ds = ms.open_output()
    out = ds.isel(lat=loc[0], lon=loc[1]).to_dataframe()[out_vars]
    good = pd.read_table('./metsim/data/three_hourly_48.3125_-120.5625',
                         names=out_vars)
    good.index = out.index

    # Make sure the data comes out right
    check_data(out, good, tol=0.2)
    ds.close()


def test_coordinate_dimension_matchup():
    """
    This test checks that MetSim correctely adds a coordinate
    if an input dataset is missing coordinate variables for the
    chunked dimensions.
    """
    var_rename = OrderedDict(
        latitude='lat', longitude='lon', mask='mask',
        elevation='elev', pptrate='prec', maxtemp='t_max', mintemp='t_min')
    filename = './examples/example_dimtest.conf'
    conf = io.read_config(DummyOpts(filename))
    conf['out_dir'] = tempfile.mkdtemp('results')
    ms = MetSim(conf)
    ds = xr.open_dataset('./metsim/data/dim_test.nc')
    assert 'hru' not in ds.coords
    assert 'hru' in ms.met_data.coords


@pytest.mark.parametrize('kind', ['ascii', 'bin', 'nc',
                                  'constant_vars_ascii',
                                  'constant_vars_bin',
                                  'constant_vars_nc'])
def test_examples(kind):
    filename = './examples/example_{kind}.conf'.format(kind=kind)
    conf = io.read_config(DummyOpts(filename))
    out_dir = tempfile.mkdtemp('results')
    conf['out_dir'] = out_dir
    ms = MetSim(conf)
    ms.run()
    assert ms.open_output() is not None


def test_yaml_config():
    filename = './examples/example_yaml.yaml'
    conf = io.read_config(DummyOpts(filename))
    out_dir = tempfile.mkdtemp('results')
    conf['out_dir'] = out_dir
    ms = MetSim(conf)
    ms.run()
    assert ms.open_output() is not None


