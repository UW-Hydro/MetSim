"""
The main object of the MetSim package. The MetSim object
is used to set up and launch forcing generation and/or
disaggregation routines.

The MetSim object uses a class dictionary to refer to
the model setup, which can be modified after instantiation
if necessary.  Before calling `run` or `launch` on the
instance it is required to call the `load` function to
ensure that all of the required parameters have been
set and that the input data is sufficient to provide the
output specified.
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import json
import logging
import itertools
import time as tm
from getpass import getuser
from multiprocessing import Pool

import numpy as np
import pandas as pd
import xarray as xr
from collections import OrderedDict, Iterable

from metsim import io
from metsim.methods import mtclim
from metsim.disaggregate import disaggregate
from metsim.physics import solar_geom
import metsim.constants as cnst
from metsim.datetime import date_range

references = '''Thornton, P.E., and S.W. Running, 1999. An improved algorithm for estimating incident daily solar radiation from measurements of temperature, humidity, and precipitation. Agricultural and Forest Meteorology, 93:211-228.
Kimball, J.S., S.W. Running, and R. Nemani, 1997. An improved method for estimating surface humidity from daily minimum temperature. Agricultural and Forest Meteorology, 85:87-98.
Bohn, T. J., B. Livneh, J. W. Oyler, S. W. Running, B. Nijssen, and D. P. Lettenmaier, 2013a: Global evaluation of MT-CLIM and related algorithms for forcing of ecological and hydrological models, Agr. Forest. Meteorol., 176, 38-49, doi:10.1016/j.agrformet.2013.03.003.'''

now = tm.ctime(tm.time())
user = getuser()
formatter = logging.Formatter(' - '.join(
    ["%(asctime)s", "%(name)s", "%(levelname)s", "%(message)s"]))

ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(formatter)
logger = logging.getLogger("metsim")

attrs = {'pet': {'units': 'mm timestep-1', 'long_name': 'potential evaporation',
                 'standard_name': 'water_potential_evaporation_flux'},
         'prec': {'units': 'mm timestep-1', 'long_name': 'precipitation',
                  'standard_name': 'precipitation_flux'},
         'shortwave': {'units': 'W m-2', 'long_name': 'shortwave radiation',
                       'standard_name': 'surface_downwelling_shortwave_flux'},
         'longwave': {'units': 'W m-2', 'long_name': 'longwave radiation',
                      'standard_name': 'surface_downwelling_longwave_flux'},
         't_max': {'units': 'C', 'long_name': 'maximum daily air temperature',
                   'standard_name': 'daily_maximum_air_temperature'},
         't_min': {'units': 'C', 'long_name': 'minimum daily air temperature',
                   'standard_name': 'daily_minimum_air_temperature'},
         'temp': {'units': 'C', 'long_name': 'air temperature',
                  'standard_name': 'air_temperature'},
         'vapor_pressure': {'units': 'kPa', 'long_name': 'vapor pressure',
                            'standard_name': 'vapor_pressure'},
         'air_pressure': {'units': 'kPa', 'long_name': 'air pressure',
                          'standard_name': 'air_pressure'},
         'tskc': {'units': 'fraction', 'long_name': 'cloud fraction',
                  'standard_name': 'cloud_fraction'},
         'rel_humid': {'units': '%', 'long_name': 'relative humidity',
                       'standard_name': 'relative_humidity'},
         'spec_humid': {'units': '', 'long_name': 'specific humidity',
                        'standard_name': 'specific_humidity'},
         'time': {'long_name': 'local time at grid location',
                  'standard_name': 'local_time'},
         '_global': {'conventions': '1.6', 'title': 'Output from MetSim',
                     'institution': 'University of Washington',
                     'source': 'metsim.py',
                     'history': 'Created: {0} by {1}'.format(now, user),
                     'references': references,
                     'documentation': 'Times given are local to their locations.',
                     'comment': 'no comment at this time'}}

attrs = {k: OrderedDict(v) for k, v in attrs.items()}


class MetSim(object):
    """
    MetSim handles the distribution of jobs that write to a common file
    by launching muliple processes and queueing up their writeback so that
    work can be done while IO is happening.
    """

    # Class variables
    methods = {'mtclim': mtclim}
    params = {
        "method": '',
        "domain": '',
        "state": '',
        "out_dir": '',
        "out_state": '',
        "out_prefix": 'forcing',
        "start": '',
        "stop": '',
        "time_step": '',
        "calendar": 'standard',
        "out_fmt": '',
        "out_precision": 'f8',
        "verbose": 0,
        "sw_prec_thresh": 0.0,
        "mtclim_swe_corr": False,
        "annual": False,
        "lw_cloud": 'cloud_deardorff',
        "lw_type": 'prata',
        "tdew_tol": 1e-6,
        "tmax_daylength_fraction": 0.67,
        "snow_crit_temp": -6.0,
        "snow_melt_rate": 0.042,
        "rain_scalar": 0.75,
        "tday_coef": 0.45,
        "lapse_rate": 0.0065,
        "iter_dims": ['lat', 'lon'],
        "out_vars": ['temp', 'prec', 'shortwave', 'longwave',
                     'vapor_pressure', 'rel_humid']
    }

    def __init__(self, params: dict):
        """
        Constructor
        """
        # Record parameters
        self.params.update(params)
        logger.setLevel(self.params['verbose'])
        ch.setLevel(self.params['verbose'])
        logger.addHandler(ch)
        self.domain = io.read_domain(self.params)
        self.met_data = io.read_met_data(self.params, self.domain)
        self.state = io.read_state(self.params, self.domain)
        self.met_data['elev'] = self.domain['elev']
        self.met_data['lat'] = self.domain['lat']
        self._aggregate_state()

    def launch(self):
        """Farm out the jobs to separate processes"""
        # Do the forcing generation and disaggregation if required
        self._validate_setup()
        nprocs = self.params['nprocs']
        time_dim = pd.DatetimeIndex(self.met_data.time.values)
        iter_list = [self.met_data[dim].values
                     for dim in self.params['iter_dims']]
        self.disagg = int(self.params['time_step']) < cnst.MIN_PER_DAY
        self.method = MetSim.methods[self.params['method']]

        if self.params['annual']:
            groups = time_dim.groupby(time_dim.year)
        else:
            groups = {'total': time_dim}

        for label, times in groups.items():
            self.pool = Pool(processes=nprocs)
            self.site_generator = itertools.product(*iter_list)
            logger.info("Beginning {}".format(label))
            status = []

            # Add in some end point data for continuity in chunking
            times_ext = times.union([times[0] - pd.Timedelta("1 days"),
                                     times[-1] + pd.Timedelta("1 days")]
                                    ).intersection(time_dim)
            data = self.met_data.sel(time=times_ext)
            self.setup_output(self.met_data.sel(time=times))
            for site in self.site_generator:
                locs = {k: v for k, v in zip(self.params['iter_dims'], site)}
                # Don't run masked cells
                if not self.domain['mask'].sel(**locs).values > 0:
                    continue
                stat = self.pool.apply_async(
                    wrap_run,
                    args=(self.method.run, locs, self.params,
                          data.sel(**locs),
                          self.state.sel(**locs),
                          self.disagg, times, label),
                    callback=self._unpack_results)
                status.append(stat)
            self.pool.close()
            # Check that everything worked
            [stat.get() for stat in status]
            self.pool.join()
            self.write(label)

    def run(self):
        """
        Kicks off the disaggregation and queues up data for IO
        """
        self._validate_setup()
        time_dim = pd.DatetimeIndex(self.met_data.time.values)
        iter_list = [self.met_data[dim].values
                     for dim in self.params['iter_dims']]
        self.disagg = int(self.params['time_step']) < cnst.MIN_PER_DAY
        self.method = MetSim.methods[self.params['method']]

        if self.params['annual']:
            groups = time_dim.groupby(time_dim.year)
        else:
            groups = {'total': time_dim}

        for label, times in groups.items():
            logger.info("Beginning {}".format(label))
            self.site_generator = itertools.product(*iter_list)
            # Add in some end point data for continuity
            times_ext = times.union([times[0] - pd.Timedelta("1 days"),
                                     times[-1] + pd.Timedelta("1 days")]
                                    ).intersection(time_dim)
            data = self.met_data.sel(time=times_ext)
            self.setup_output(self.met_data.sel(time=times))

            for site in self.site_generator:
                locs = {k: v for k, v in zip(self.params['iter_dims'], site)}
                # Don't run masked cells
                if not self.domain['mask'].sel(**locs).values > 0:
                    continue
                wrap_results = wrap_run(
                    self.method.run, locs, self.params, data.sel(**locs),
                    self.state.sel(**locs), self.disagg, times, label)

                # Cut the returned data down to the correct time index
                self._unpack_results(wrap_results)

            self.write(label)

    def _unpack_results(self, result: tuple):
        """Put results into the master dataset"""
        if len(result) == 3:
            locs, df, state = result
            self._unpack_state(state, locs)
        else:
            locs, df = result
        for varname in self.params['out_vars']:
            try:
                self.output[varname].loc[locs] = df[varname]
            except ValueError as e:
                logger.error(e)
                logger.error("This error is probably indicitive of a mismatch "
                             "between the domain and input data. Check that "
                             "all of your cells inside of the mask have both "
                             "elevation in the domain as well as all of the "
                             "required input forcings.")
                raise

    def _unpack_state(self, result: pd.DataFrame, locs: dict):
        """Put restart values in the state dataset"""
        # We concatenate with the old state values in case we don't
        # have 90 new days to use
        tmin = np.concatenate((self.state['t_min'].sel(**locs).values[:],
                               result['t_min'].values))
        tmax = np.concatenate((self.state['t_max'].sel(**locs).values[:],
                               result['t_max'].values))
        prec = np.concatenate((self.state['prec'].sel(**locs).values[:],
                               result['prec'].values))
        self.state['t_min'].sel(**locs).values[:] = tmin[-90:]
        self.state['t_max'].sel(**locs).values[:] = tmax[-90:]
        self.state['prec'].sel(**locs).values[:] = prec[-90:]
        self.state['swe'].sel(**locs).values = result['swe'].values[-1]
        state_start = result.index[-1] - pd.Timedelta('89 days')
        self.state.time.values = date_range(state_start, result.index[-1],
                                            calendar=self.params['calendar'])

    def setup_output(self, prototype: xr.Dataset=None):
        if not prototype:
            prototype = self.met_data
        self.disagg = int(self.params['time_step']) < cnst.MIN_PER_DAY
        # Number of timesteps
        if self.disagg:
            delta = pd.Timedelta('1 days') - pd.Timedelta(
                '{} minutes'.format(self.params['time_step']))
        else:
            delta = pd.Timedelta('0 days')

        start = pd.Timestamp(prototype.time.values[0]).to_pydatetime()
        stop = pd.Timestamp(prototype.time.values[-1]).to_pydatetime()
        times = date_range(start, stop + delta,
                           freq="{}T".format(self.params['time_step']),
                           calendar=self.params['calendar'])
        n_ts = len(times)

        shape = (n_ts, ) + self.domain['mask'].shape
        dims = ('time', ) + self.domain['mask'].dims
        coords = {'time': times, **self.domain.mask.coords}
        self.output = xr.Dataset(coords=coords)
        if 'elev' in self.params:
            self.params.pop('elev')
        for k, v in self.params.items():
            # Need to convert some parameters to strings
            if k in ['start', 'stop', 'annual', 'mtclim_swe_corr']:
                v = str(v)
            # Don't include complex types
            if isinstance(v, dict):
                v = json.dumps(v)
            elif not isinstance(v, str) and isinstance(v, Iterable):
                v = ', '.join(v)
            attrs['_global'][k] = v
        self.output.attrs = attrs['_global']
        for varname in self.params['out_vars']:
            self.output[varname] = xr.DataArray(
                data=np.full(shape, np.nan),
                coords=coords, dims=dims,
                name=varname, attrs=attrs.get(varname, {}),
                encoding={'dtype': self.params['out_precision'],
                          '_FillValue': cnst.FILL_VALUES['f8']})
        self.output.time.attrs.update(attrs['time'])

    def _aggregate_state(self):
        """Aggregate data out of the state file and load it into `met_data`"""
        # Precipitation record
        begin_record = self.params['start'] - pd.Timedelta("90 days")
        end_record = self.params['start'] - pd.Timedelta("1 days")
        record_dates = date_range(begin_record, end_record,
                                  calendar=self.params['calendar'])
        trailing = self.state['prec'].sel(time=record_dates)
        total_precip = xr.concat([trailing, self.met_data['prec']], dim='time')
        total_precip = (cnst.DAYS_PER_YEAR * total_precip.rolling(
            time=90).mean().drop(record_dates, dim='time'))
        self.met_data['seasonal_prec'] = total_precip

        # Smoothed daily temperature range
        trailing = self.state['t_max'] - self.state['t_min']
        begin_record = self.params['start'] - pd.Timedelta("90 days")
        end_record = self.params['start'] - pd.Timedelta("1 days")
        record_dates = date_range(begin_record, end_record,
                                  calendar=self.params['calendar'])
        trailing['time'] = record_dates
        dtr = self.met_data['t_max'] - self.met_data['t_min']
        sm_dtr = xr.concat([trailing, dtr], dim='time')
        sm_dtr = sm_dtr.rolling(time=30).mean().drop(record_dates, dim='time')
        self.met_data['smoothed_dtr'] = sm_dtr

        # Put in SWE data
        self.state['swe'] = self.state.sel(time=end_record).swe.drop('time')
        self.met_data['swe'] = self.state.swe.copy()

    def _validate_setup(self):
        """Updates the global parameters dictionary"""
        errs = [""]

        # Make sure there's some input
        if not len(self.params['forcing']):
            errs.append("Requires input forcings to be specified")

        # Parameters that can't be empty strings or None
        non_empty = ['method', 'out_dir', 'out_state', 'start',
                     'stop', 'time_step', 'out_fmt',
                     'forcing_fmt', 'domain_fmt', 'state_fmt']
        for each in non_empty:
            if self.params[each] is None or self.params[each] == '':
                errs.append("Cannot have empty value for {}".format(each))

        # Make sure time step divides evenly into a day
        if cnst.MIN_PER_DAY % int(self.params['time_step']):
            errs.append("Time step must divide 1440 evenly.  Got {}"
                        .format(self.params['time_step']))

        # Check for required input variable specification
        required_in = ['t_min', 't_max', 'prec']
        for each in required_in:
            if each not in self.met_data.variables:
                errs.append("Input requires {}".format(each))

        # Make sure that we are going to write out some data
        if not len(self.params['out_vars']):
            errs.append("Output variable list must not be empty")

        # Check that the parameters specified are available
        opts = {'mtclim_swe_corr': [True, False],
                'out_precision': ['f4', 'f8'],
                'lw_cloud': ['default', 'cloud_deardorff'],
                'lw_type': ['default', 'tva', 'anderson',
                            'brutsaert', 'satterlund',
                            'idso', 'prata']}
        for k, v in opts.items():
            if not self.params[k] in v:
                errs.append("Invalid option given for {}".format(k))

        # If any errors, raise and give a summary
        if len(errs) > 1:
            raise Exception("\n  ".join(errs))

    def write(self, suffix=""):
        """
        Dispatch to the right function based on the configuration given
        """
        dispatch = {
                'netcdf': self.write_netcdf,
                'ascii': self.write_ascii,
                'data': self.write_data
                }
        dispatch[self.params.get('out_fmt', 'netcdf').lower()](suffix)

    def write_netcdf(self, suffix: str):
        """Write out as NetCDF to the output file"""
        logger.info("Writing netcdf...")
        for dirname in [self.params['out_dir'],
                        os.path.dirname(self.params['out_state'])]:
            os.makedirs(dirname, exist_ok=True)

        # all state variables are written as doubles
        state_encoding = {'time': {'dtype': 'f8'}}
        for v in self.state:
            state_encoding[v] = {'dtype': 'f8'}
        # write state file

        self.state.to_netcdf(self.params['out_state'], encoding=state_encoding)
        # write output file
        fname = '{}_{}.nc'.format(self.params['out_prefix'], suffix)
        output_filename = os.path.join(self.params['out_dir'], fname)
        self.output.to_netcdf(output_filename,
                              unlimited_dims=['time'],
                              encoding={'time': {'dtype': 'f8'}})

    def write_ascii(self, suffix):
        """Write out as ASCII to the output file"""
        logger.info("Writing ascii...")
        for dirname in [self.params['out_dir'],
                        os.path.dirname(self.params['out_state'])]:
            os.makedirs(dirname, exist_ok=True)
        # all state variables are written as doubles
        state_encoding = {'time': {'dtype': 'f8'}}
        for v in self.state:
            state_encoding[v] = {'dtype': 'f8'}
        # write state file
        self.state.to_netcdf(self.params['out_state'], encoding=state_encoding)
        # Need to create new generator to loop over
        iter_list = [self.met_data[dim].values
                     for dim in self.params['iter_dims']]
        site_generator = itertools.product(*iter_list)

        for site in site_generator:
            locs = {k: v for k, v in zip(self.params['iter_dims'], site)}
            if not self.domain['mask'].sel(**locs).values > 0:
                continue
            fname = ("{}_" * (len(iter_list)+1) + "{}.csv").format(
                self.params['out_prefix'], *site, suffix)
            fpath = os.path.join(self.params['out_dir'], fname)
            self.output.sel(**locs)[self.params[
                   'out_vars']].to_dataframe().to_csv(fpath)

    def write_data(self, suffix):
        pass


def wrap_run(func: callable, loc: dict, params: dict,
             ds: xr.Dataset, state: xr.Dataset, disagg: bool,
             out_times: pd.DatetimeIndex, year: str):
    """
    Iterate over a chunk of the domain. This is wrapped
    so we can return a tuple of locs and df.

    Parameters
    ----------
    func: callable
        The function to call to do the work
    loc: dict
        Some subset of the domain to do work on
    params: dict
        Parameters from a MetSim object
    ds: xr.Dataset
        Input forcings and domain
    state: xr.Dataset
        State variables at the point of interest
    disagg: bool
        Whether or not we should run a disagg routine
    out_times: pd.DatetimeIndex
        Times to return (should be trimmed 1 day at
        each end from the given index)
    year: str
        The year being run. This is used to add on
        extra times to make output smooth at endpoints
        if the run is chunked in time.

    Returns
    -------
    results
        A list of tuples arranged as
        (location, hourly_output, daily_output)
    """
    logger.info("Processing {}".format(loc))
    lat = ds['lat'].values
    elev = ds['elev'].values
    swe = ds['swe'].values
    params['elev'] = elev
    df = ds.to_dataframe()
    # solar_geom returns a tuple due to restrictions of numba
    # for clarity we convert it to a dictionary here
    sg = solar_geom(elev, lat, params['lapse_rate'])
    sg = {'tiny_rad_fract': sg[0], 'daylength': sg[1],
          'potrad': sg[2], 'tt_max0': sg[3]}

    # Generate the daily values - these are saved
    # so that we can use a subset of them to write
    # out the state file later
    df_base = func(df, params, sg, elev=elev, swe=swe)

    if disagg:
        # Get some values for padding the time list,
        # so that when interpolating in the disaggregation
        # functions we can match endpoints with adjoining
        # chunks - if no data is available, just repeat some
        # default values (this case is used at the very
        # beginning and end of the record)
        try:
            prevday = out_times[0] - pd.Timedelta('1 days')
            t_begin = [ds['t_min'].sel(time=prevday),
                       ds['t_max'].sel(time=prevday)]
        except (KeyError, ValueError):
            t_begin = [state['t_min'].values[-1],
                       state['t_max'].values[-1]]
        try:
            nextday = pd.datetime(int(year)+1, 1, 1)
            t_end = [ds['t_min'].sel(time=nextday),
                     ds['t_max'].sel(time=nextday)]
        except (KeyError, ValueError):
            # None so that we don't extend the record
            t_end = None

        # Disaggregate to subdaily values
        df_complete = disaggregate(df, params, sg, t_begin, t_end)
        # Calculate the times that we want to get out by chopping
        # off the endpoints that were added on previously
        start = out_times[0]
        stop = (out_times[-1] + pd.Timedelta('1 days')
                - pd.Timedelta("{} minutes".format(params['time_step'])))
        new_times = date_range(
            start, stop, freq='{}T'.format(params['time_step']),
            calendar=params['calendar'])
    else:
        # convert srad to daily average flux from daytime flux
        df_base['shortwave'] *= df_base['dayl'] / cnst.SEC_PER_DAY
        # If we're outputting daily values, we dont' need to
        # change the output dates - see inside of `if` condition
        # above for more explanation
        new_times = out_times
        df_complete = df_base

    # Cut the returned data down to the correct time index
    df_complete = df_complete.loc[new_times[0]:new_times[-1]]
    return (loc, df_complete, df_base)
