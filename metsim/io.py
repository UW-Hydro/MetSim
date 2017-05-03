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

import numpy as np
import pandas as pd
import xarray as xr


def read_met_data(params, domain):
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


def process_nc(params, domain):
    read_funcs = {
        "netcdf": read_netcdf,
        "data": read_data
    }
    return read_funcs[params['forcing_fmt']](
        params['forcing'], params['start'],
        params['stop'], params.get('forcing_vars', None))


def process_vic(params, domain):
    """Process all files to find spatial extent"""
    read_funcs = {
        "binary": read_binary,
        "ascii": read_ascii,
    }
    # Creates the master dataset which will be used to parallelize
    met_data = xr.Dataset(
        coords={'time': params['dates'],
                'lon': domain['lon'],
                'lat': domain['lat']},
        attrs={'n_days': len(params['dates'])})
    shape = (len(params['dates']),
             len(domain['lat']),
             len(domain['lon']))

    for var in params['forcing_vars']:
        met_data[var] = xr.Variable(
            ('time', 'lat', 'lon'), np.full(shape, np.nan))

    # Fill in the data
    ilist, jlist = np.nonzero(np.nan_to_num(domain['mask'].values))
    for job in params['forcing']:
        _, lat, lon = os.path.basename(job).split("_")
        i = np.unique(ilist)[
            list(domain['lat'].values).index(float(lat))]
        j = np.unique(jlist)[
            list(domain['lon'].values).index(float(lon))]
        ds = read_funcs[params['forcing_fmt']](
            job, params['start'], params['stop'], params['forcing_vars'])
        for var in params['forcing_vars'].keys():
            met_data[var].values[:, i, j] = ds[var].values


def read_domain(params):
    read_funcs = {
        "netcdf": read_netcdf,
        "data": read_data
    }
    return read_funcs[params['domain_fmt']](
        params['domain'], None, None, params.get('domain_vars', None))


def read_state(params):
    read_funcs = {
        "netcdf": read_netcdf,
        "data": read_data
    }
    start = params['start'] - pd.Timedelta('90 days')
    stop = params['start'] - pd.Timedelta('1 days')
    return read_funcs[params['state_fmt']](
        params['state'], start, stop, params.get('state_vars', None))


def read_ascii(data_handle, start, stop, var_dict) -> xr.Dataset:
    """Read in an ascii forcing file"""
    dates = pd.date_range(start, stop)
    names = var_dict.keys()
    ds = pd.read_table(data_handle, header=None, delim_whitespace=True,
                       names=names).head(len(dates))
    ds.index = dates
    return ds


def read_netcdf(data_handle, start, stop, var_dict) -> xr.Dataset:
    """
    Read in a NetCDF file and add elevation information
    """
    ds = xr.open_dataset(data_handle)
    if var_dict is not None:
        ds.rename(var_dict, inplace=True)
        varlist = list(var_dict.values())

    if start is not None and stop is not None:
        ds = ds[varlist].sel(time=slice(start, stop))
        dates = ds.indexes['time']
        ds['day_of_year'] = xr.Variable(('time', ), dates.dayofyear)

    return ds.load()


def read_data(data_handle, start, stop, var_dict):
    varlist = list(data_handle.keys())
    if var_dict is not None:
        data_handle.rename(var_dict, inplace=True)
        varlist = list(var_dict.values())

    if start is not None and stop is not None:
        data_handle = data_handle[varlist].sel(time=slice(start, stop))
        dates = data_handle.indexes['time']
        data_handle['day_of_year'] = xr.Variable(('time', ), dates.dayofyear)

    return data_handle.load()


def read_binary(data_handle, start, stop, var_dict) -> xr.Dataset:
    """ Reads a binary forcing file (VIC 4 format) """
    dates = pd.date_range(start, stop)
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
        points_needed = 4*n_days
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
    param_list = os.path.basename(data_handle).split("_")
    params = {"name": param_list[0],
              "lat": float(param_list[1]),
              "lon": float(param_list[2]),
              "n_days": int(n_days)}

    # Assemble the dataset
    data_dict = {c[0]: (['time'], c[1]) for c in zip(var_names, data_list)}
    data_dict['elev'] = (['lon', 'lat'], params['elev'])
    data_dict['day_of_year'] = (['time'], dates.dayofyear)
    df = xr.Dataset(data_dict,
                    coords={'lon': [params['lon']],
                            'lat': [params['lat']],
                            'time': dates},
                    attrs={'n_days': params['n_days']})
    return df
