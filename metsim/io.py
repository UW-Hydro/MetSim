"""
IO Module for MetSim
"""
# Meteorology Simulator
# Copyright (C) 2017  The Computational Hydrology Group, Department of Civil
# and Environmental Engineering, University of Washington.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

import os
import struct
import logging
import json
import yaml
import warnings

import numpy as np
import pandas as pd
import xarray as xr

from collections import OrderedDict
from configparser import ConfigParser
from metsim.datetime import date_range
from metsim import metsim


def _invert_dict(d):
    return OrderedDict([reversed(item) for item in d.items()])


def _to_list(s):
    return json.loads(s.replace("'", '"').split('#')[0])


def check_config_type(config_file):
    # Search for first non-comment line
    line = '#'
    with open(config_file, 'r') as f:
        while line.strip().startswith('#'):
            line = f.readline()
    # Strip out any inline comment
    line = line.split('#')[0].strip()
    if line == '[MetSim]':
        return 'ini'
    if line == 'MetSim:':
        return 'yaml'
    # fallback
    return 'ini'


def read_config(opts):
    """Initialize some information based on the options & config"""
    config = ConfigParser()
    config.optionxform = str
    config_file = opts.config
    config_type = check_config_type(config_file)
    config_readers = {'yaml': read_yaml_config,
                      'ini': read_ini_config}
    return config_readers[config_type](config_file, opts)


def read_yaml_config(config_file, opts):
    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
    conf = OrderedDict(config['MetSim'])
    conf['forcing_vars'] = _invert_dict(OrderedDict(config['forcing_vars']))
    conf['domain_vars'] = _invert_dict(OrderedDict(config['domain_vars']))
    conf['state_vars'] = _invert_dict(OrderedDict(config['state_vars']))
    conf['out_vars'] = OrderedDict(config['out_vars'])
    for k, v in conf['out_vars'].items():
        if 'units' not in v:
            v['units'] = metsim.available_outputs[k]['units']
    conf['chunks'] = OrderedDict(config['chunks'])
    if 'constant_vars' in config:
        conf['constant_vars'] = OrderedDict(config['constant_vars'])

    # If the forcing variable is a directory, scan it for files
    if os.path.isdir(conf['forcing']):
        forcing_files = [os.path.join(conf['forcing'], fn) for fn in
                         next(os.walk(conf['forcing']))[2]]
    else:
        forcing_files = conf['forcing']

    # Update the full configuration
    conf.update({"calendar": conf.get('calendar', 'standard'),
                 "scheduler": opts.scheduler,
                 "num_workers": opts.num_workers,
                 "verbose": logging.DEBUG if opts.verbose else logging.INFO,
                 "forcing": forcing_files,
                 "out_dir": os.path.abspath(conf['out_dir'])})

    conf = {k: v for k, v in conf.items() if v != []}
    return conf


def read_ini_config(config_file, opts):
    config = ConfigParser()
    config.optionxform = str
    config.read(config_file)
    conf = OrderedDict(config['MetSim'])
    conf['forcing_vars'] = OrderedDict(config['forcing_vars'])
    if conf['forcing_fmt'] != 'binary':
        conf['forcing_vars'] = _invert_dict(conf['forcing_vars'])
    conf['domain_vars'] = _invert_dict(OrderedDict(config['domain_vars']))
    conf['state_vars'] = _invert_dict(OrderedDict(config['state_vars']))
    conf['chunks'] = OrderedDict(config['chunks'])
    if 'constant_vars' in config:
        conf['constant_vars'] = OrderedDict(config['constant_vars'])

    # If the forcing variable is a directory, scan it for files
    if os.path.isdir(conf['forcing']):
        forcing_files = [os.path.join(conf['forcing'], fn) for fn in
                         next(os.walk(conf['forcing']))[2]]
    else:
        forcing_files = conf['forcing']

    # Ensure that parameters with boolean values are correctly recorded
    for bool_param in ['utc_offset', 'period_ending']:
        if (bool_param in conf.keys()
            and conf[bool_param].strip().lower() == 'true'):
            conf[bool_param] = True
        else:
            conf[bool_param] = False

    # Update the full configuration
    conf.update({"calendar": conf.get('calendar', 'standard'),
                 "scheduler": opts.scheduler,
                 "num_workers": opts.num_workers,
                 "verbose": logging.DEBUG if opts.verbose else logging.INFO,
                 "forcing": forcing_files,
                 "out_dir": os.path.abspath(conf['out_dir']),
                 "prec_type": conf.get('prec_type', 'uniform')})

    # List variant
    if 'out_vars' in conf:
        conf['out_vars'] = _to_list(conf['out_vars'])
        temp = {}
        for ov in conf['out_vars']:
            temp[ov] = metsim.available_outputs[ov]
        conf['out_vars'] = temp
    else:
        conf['out_vars'] = {}
    # Dict variant
    if 'out_vars' in config:
        temp = {}
        for varname, outname in config['out_vars'].items():
            temp[varname] = metsim.available_outputs[varname]
            temp[varname]['out_name'] = outname
        conf['out_vars'].update(temp)

    conf = {k: v for k, v in conf.items() if v != []}
    return conf



def read_met_data(params: dict, domain: xr.Dataset) -> xr.Dataset:
    """
    Read input meteorological forcings for MetSim.
    This method supports ascii, binary, netcdf, and
    xarray input pointers.  The input source is derived
    from the key 'forcing' in the params dictionary.
    The format of the data is derived from 'in_format'
    key in the parameter dictionary.
    """
    process_funcs = {
        "netcdf": process_nc,
        "binary": process_vic,
        "ascii": process_vic,
        "data": process_nc
    }
    return process_funcs[params['forcing_fmt']](params, domain)


def read_domain(params: dict) -> xr.Dataset:
    """Load in a domain file"""
    return read_netcdf(
        params['domain'], calendar=params['calendar'],
        var_dict=params.get('domain_vars', None))


def read_state(params: dict, domain: xr.Dataset) -> xr.Dataset:
    """Load in a state file"""
    return read_netcdf(
        params['state'], domain=domain,
        start=params['state_start'], stop=params['state_stop'],
        calendar=params['calendar'],
        var_dict=params.get('state_vars', None))


def process_nc(params: dict, domain: xr.Dataset) -> xr.Dataset:
    """Process NetCDF-like Data"""
    read_funcs = {
        "netcdf": read_netcdf,
        "data": read_data
    }
    return read_funcs[params['forcing_fmt']](
        params['forcing'], domain=domain,
        start=params['start'], stop=params['stop'],
        calendar=params['calendar'], var_dict=params.get('forcing_vars', None))


def process_vic(params: dict, domain: xr.Dataset) -> xr.Dataset:
    """Process VIC-like data"""
    read_funcs = {
        "binary": read_binary,
        "ascii": read_ascii,
    }

    # Creates the master dataset which will be used to parallelize
    dates = date_range(params['start'], params['stop'],
                       calendar=params['calendar'])
    coords = {'time': dates,
              'lon': domain['lon'],
              'lat': domain['lat']}
    shape = (len(dates), len(domain['lat']), len(domain['lon']))
    dims = ('time', 'lat', 'lon', )

    met_data = xr.Dataset(
        coords=coords, attrs={'n_days': len(dates)})
    for var in params['forcing_vars']:
        met_data[var] = xr.DataArray(data=np.full(shape, np.nan),
                                     coords=coords, dims=dims, name=var)

    # Fill in the data
    for job in params['forcing']:
        try:
            _, lat, lon = os.path.basename(job).split("_")[-3:]
            lat, lon = float(lat), float(lon)
            if not domain['mask'].sel(lat=lat, lon=lon).values > 0:
                continue
            ds = read_funcs[params['forcing_fmt']](
                job, start=params['start'], stop=params['stop'],
                calendar=params['calendar'], var_dict=params['forcing_vars'])
            for var in params['forcing_vars'].keys():
                met_data[var].loc[{'lat': lat, 'lon': lon}] = ds[var]
        except (ValueError, KeyError):
            continue
    return met_data


def read_ascii(data_handle, domain=None, is_worker=False,
               start=None, stop=None, calendar='standard',
               var_dict=None) -> xr.Dataset:
    """Read in an ascii forcing file"""
    dates = date_range(start, stop, calendar=calendar)
    names = list(var_dict.keys())
    ds = pd.read_csv(data_handle, header=None, delim_whitespace=True,
                     names=names).head(len(dates))
    ds.index = dates
    return ds


def read_netcdf(data_handle, domain=None, is_worker=False,
                start=None, stop=None, calendar='standard',
                var_dict=None) -> xr.Dataset:
    """Read in a NetCDF file"""
    if '*' in data_handle:
        ds = xr.open_mfdataset(data_handle)
    else:
        ds = xr.open_dataset(data_handle)

    if domain is not None:
        ds = ds.sel({k: domain[k]
                     for k in list(domain.dims.keys())
                     if k in list(ds.dims.keys())})
    else:
        dims_wo_coords = set(ds.dims) - set(ds.coords)
        for d in dims_wo_coords:
            if is_worker:
                logger = logging.getLogger('MetSim')
                logger.warning(
                    'Setting sequential coordinate on dimension {}'.format(d))
            ds[d] = np.arange(0, len(ds[d]))

    if 'time' in ds.coords:
        if isinstance(ds.indexes['time'], xr.CFTimeIndex):
            ds['time'] = ds.indexes['time'].to_datetimeindex()
        ds['time'] = (ds.indexes['time'] -
                      pd.Timedelta(hours=11, minutes=59, seconds=59)).round('D')

    if var_dict is not None:
        var_list = list(var_dict.keys())
        ds = ds[var_list]
        ds = ds.rename(var_dict)

    if start is not None or stop is not None:
        ds = ds.sel(time=slice(start, stop))
        dates = ds.indexes['time']
        ds['day_of_year'] = xr.Variable(('time', ), dates.dayofyear)

    return ds


def read_data(data_handle, domain=None, is_worker=False,
              start=None, stop=None, calendar='standard',
              var_dict=None) -> xr.Dataset:
    """Read data directly from an xarray dataset"""
    varlist = list(data_handle.keys())
    if var_dict is not None:
        data_handle = data_handle.rename(var_dict)
        varlist = list(var_dict.values())
    data_handle = data_handle[varlist]

    if start is not None and stop is not None:
        data_handle = data_handle.sel(time=slice(start, stop))
        dates = data_handle.indexes['time']
        data_handle['day_of_year'] = xr.Variable(('time', ), dates.dayofyear)

    return data_handle


def read_binary(data_handle, domain=None,
                start=None, stop=None, calendar='standard',
                var_dict=None) -> xr.Dataset:
    """Reads a binary forcing file (VIC 4 format)"""
    dates = date_range(start, stop, calendar=calendar)
    n_days = len(dates)
    type_lookup = {'signed': 'h', 'unsigned': 'H'}
    # Pack these for nicer syntax in the loop
    var_names = var_dict.keys()
    data_list = [[] for var in var_dict.keys()]
    n_vars = len(var_names)
    params = [p for p in var_dict.values()]
    scales = [float(s.split()[0]) for s in params]
    datatypes = [type_lookup[s.split()[-1]] for s in params]
    with open(data_handle, 'rb') as f:
        i = 0
        points_read = 0
        points_needed = 4 * n_days
        while points_read != points_needed:
            bytes = f.read(2)
            if bytes:
                # Get correct variable and data type with i,
                # then unpack & scale
                data_list[i].append(
                    struct.unpack(datatypes[i], bytes)[0] / scales[i])
                i = (i + 1) % n_vars
                points_read += 1
            else:
                break

    # Binary forcing files have naming format $NAME_$LAT_$LON
    param_list = os.path.basename(data_handle).split("_")[-3:]
    params = {"name": param_list[0],
              "lat": float(param_list[1]),
              "lon": float(param_list[2]),
              "n_days": int(n_days)}

    # Assemble the dataset
    data_dict = {c[0]: (['time'], c[1]) for c in zip(var_names, data_list)}
    data_dict['day_of_year'] = (['time'], dates.dayofyear)
    df = xr.Dataset(data_dict,
                    coords={'lon': [params['lon']],
                            'lat': [params['lat']],
                            'time': dates},
                    attrs={'n_days': params['n_days']})
    return df
