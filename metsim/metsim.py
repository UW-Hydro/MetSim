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

import itertools
import json
import logging
import os
import sys
import time as tm
from collections import Iterable, OrderedDict
from getpass import getuser

import numpy as np
import pandas as pd
import xarray as xr

import dask
from dask.diagnostics import ProgressBar
from netCDF4 import Dataset
from cftime import date2num

from xarray.backends.locks import get_write_lock, combine_locks, NETCDFC_LOCK

import metsim.constants as cnst
from metsim import io
from metsim.datetime import date_range
from metsim.disaggregate import disaggregate
from metsim.methods import mtclim
from metsim.physics import solar_geom

NO_SLICE = {}
DASK_CORE_SCHEDULERS = ['multiprocessing', 'threading', 'synchronous',
                        'processes', 'threads', 'single-threaded', 'sync']

references = '''Thornton, P.E., and S.W. Running, 1999. An improved algorithm for estimating incident daily solar radiation from measurements of temperature, humidity, and precipitation. Agricultural and Forest Meteorology, 93:211-228.
Kimball, J.S., S.W. Running, and R. Nemani, 1997. An improved method for estimating surface humidity from daily minimum temperature. Agricultural and Forest Meteorology, 85:87-98.
Bohn, T. J., B. Livneh, J. W. Oyler, S. W. Running, B. Nijssen, and D. P. Lettenmaier, 2013a: Global evaluation of MT-CLIM and related algorithms for forcing of ecological and hydrological models, Agr. Forest. Meteorol., 176, 38-49, doi:10.1016/j.agrformet.2013.03.003.'''

DEFAULT_TIME_UNITS = 'hours since 2000-01-01 00:00:00.0'

now = tm.ctime(tm.time())
user = getuser()

attrs = {'pet': {'units': 'mm timestep-1',
                 'long_name': 'potential evaporation',
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
         'wind': {'units': 'm/s', 'long_name': 'wind speed',
                  'standard_name': 'wind_speed'},
         '_global': {'conventions': '1.6', 'title': 'Output from MetSim',
                     'institution': 'University of Washington',
                     'source': 'metsim.py',
                     'history': 'Created: {0} by {1}'.format(now, user),
                     'references': references,
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
        "method": 'mtclim',
        "domain": '',
        "state": '',
        "out_dir": '',
        "out_prefix": 'forcing',
        "start": 'forcing',
        "stop": 'forcing',
        "time_step": -1,
        "calendar": 'standard',
        "prec_type": 'uniform',
        "out_precision": 'f4',
        "verbose": 0,
        "sw_prec_thresh": 0.0,
        "utc_offset": False,
        "lw_cloud": 'cloud_deardorff',
        "lw_type": 'prata',
        "tdew_tol": 1e-6,
        "tmax_daylength_fraction": 0.67,
        "rain_scalar": 0.75,
        "tday_coef": 0.45,
        "lapse_rate": 0.0065,
        "out_vars": ['temp', 'prec', 'shortwave', 'longwave',
                     'vapor_pressure', 'rel_humid'],
        "out_freq": None,
        "chunks": NO_SLICE,
        "scheduler": 'distributed',
        "num_workers": 1,
    }

    def __init__(self, params: dict, domain_slice=NO_SLICE):
        """
        Constructor
        """
        self._domain = None
        self._met_data = None
        self._state = None
        self._client = None
        self._domain_slice = domain_slice
        self.progress_bar = ProgressBar()
        self.params.update(params)
        logging.captureWarnings(True)
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(self.params['verbose'])

        formatter = logging.Formatter(' - '.join(
            ['%asctime)s', '%(name)s', '%(levelname)s', '%(message)s']))
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(formatter)
        ch.setLevel(self.params['verbose'])
        # set global dask scheduler
        if domain_slice is NO_SLICE:
            if self.params['scheduler'] in DASK_CORE_SCHEDULERS:
                dask.config.set(scheduler=self.params['scheduler'])
            else:
                from distributed import Client, progress
                if 'distributed' == self.params['scheduler']:
                    self._client = Client(
                        n_workers=self.params['num_workers'],
                        threads_per_worker=1)
                    if self.params['verbose'] == logging.DEBUG:
                        self.progress_bar = progress
                elif os.path.isfile(self.params['scheduler']):
                    self._client = Client(
                        scheduler_file=self.params['scheduler'])
                else:
                    self._client = Client(self.params['scheduler'])
        else:
            dask.config.set(scheduler=self.params['scheduler'])

        # Set up logging
        # If in verbose mode set up the progress bar
        if self.params['verbose'] == logging.DEBUG:
            if 'distributed' != self.params['scheduler']:
                self.progress_bar.register()
                self.progress_bar = lambda x: x
        else:
            # If not in verbose mode, create a dummy function
            self.progress_bar = lambda x: x
        # Create time vector(s)
        self._times = self._get_output_times(freq=self.params['out_freq'])

    def _validate_force_times(self, force_times):

        for p, i in [('start', 0), ('stop', -1)]:
            # infer times from force_times
            if isinstance(self.params[p], str):
                if self.params[p] == 'forcing':
                    self.params[p] = pd.Timestamp(
                        force_times.values[i]).to_pydatetime()
                elif '/' in self.params[p]:
                    year, month, day = map(int, self.params[p].split('/'))
                    self.params[p] = pd.datetime(year, month, day)
                else:
                    self.params[p] = pd.to_datetime(self.params[p])

        # update calendar from input data (fall back to params version)
        self.params['calendar'] = self.met_data['time'].encoding.get(
            'calendar', self.params['calendar'])

        assert self.params['start'] >= pd.Timestamp(
            force_times.values[0]).to_pydatetime()
        assert self.params['stop'] <= pd.Timestamp(
            force_times.values[-1]).to_pydatetime()

        self.params['state_start'] = (self.params['start'] -
                                      pd.Timedelta("90 days"))
        self.params['state_stop'] = (self.params['start'] -
                                     pd.Timedelta("1 days"))
        if self.params['utc_offset']:
            attrs['time'] = {'units': DEFAULT_TIME_UNITS,
                             'long_name': 'UTC time',
                             'standard_name': 'utc_time'}
        else:
            attrs['time'] = {'units': DEFAULT_TIME_UNITS,
                             'long_name': 'local time at grid location',
                             'standard_name': 'local_time'}

    @property
    def domain(self):
        if self._domain is None:
            self._domain = io.read_domain(self.params).isel(
                **self._domain_slice)
        return self._domain

    @property
    def met_data(self):
        if self._met_data is None:
            if self.domain is None:
                self._domain = io.read_domain(self.params).isel(
                    **self._domain_slice)
            self._met_data = io.read_met_data(self.params, self._domain)
            self._met_data['elev'] = self.domain['elev']
            self._met_data['lat'] = self.domain['lat']
            self._met_data['lon'] = self.domain['lon']

            # process constant_vars
            constant_vars = self.params.get('constant_vars', None)
            if constant_vars:
                da_template = self._met_data[list(self._met_data)[0]]
                for var in constant_vars.keys():
                    self._met_data[var] = xr.full_like(da_template,
                                                       float(constant_vars[var]))

            self._validate_force_times(force_times=self._met_data['time'])
        return self._met_data

    @property
    def state(self):
        if self._state is None:
            self._state = io.read_state(self.params, self.domain)
            self._aggregate_state()
        return self._state

    @property
    def slices(self):
        if not self.params['chunks']:
            return [{d: slice(None) for d in self.domain[['mask']].dims}]
        return chunk_domain(self.params['chunks'], self.domain[['mask']].dims)

    def open_output(self):
        filenames = [self._get_output_filename(times) for times in self._times]
        return xr.open_mfdataset(filenames)

    def run(self):
        self._validate_setup()
        write_locks = {}
        for times in self._times:
            filename = self._get_output_filename(times)
            self.setup_netcdf_output(filename, times)
            write_locks[filename] = combine_locks([NETCDFC_LOCK, get_write_lock(filename)])
        self.logger.info('Starting {} chunks...'.format(len(self.slices)))

        delayed_objs = [wrap_run_slice(self.params, write_locks, dslice)
                        for dslice in self.slices]
        persisted = dask.persist(delayed_objs, num_workers=self.params['num_workers'])
        self.progress_bar(persisted)
        dask.compute(persisted)
        self.logger.info('Cleaning up...')
        try:
            self._client.cluster.close()
            self._client.close()
            if self.params['verbose'] == logging.DEBUG:
                print('closed dask cluster/client')
        except Exception:
            pass

    def load_inputs(self, close=True):
        self._domain = self.domain.load()
        self._met_data = self.met_data.load()
        self._state = self.state.load()
        if close:
            self._domain.close()
            self._met_data.close()
            self._state.close()

    def setup_netcdf_output(self, filename, times):
        '''setup a single netcdf file'''
        with Dataset(filename, mode="w") as ncout:
            # dims
            dim_sizes = (None, ) + self.domain['mask'].shape
            var_dims = ('time', ) + self.domain['mask'].dims
            chunksizes = [len(times)]
            for d, s in zip(var_dims[1:], dim_sizes[1:]):
                c = int(self.params['chunks'].get(d, s))
                if c <= s:
                    chunksizes.append(c)
                else:
                    chunksizes.append(s)
            create_kwargs = {'chunksizes': chunksizes}
            for d, size in zip(var_dims, dim_sizes):
                ncout.createDimension(d, size)
            # vars
            for varname in self.params['out_vars']:
                ncout.createVariable(
                    varname, self.params['out_precision'], var_dims,
                    **create_kwargs)

            # add metadata and coordinate variables (time/lat/lon)
            time_var = ncout.createVariable('time', 'i4', ('time', ))
            time_var.calendar = self.params['calendar']
            time_var[:] = date2num(
                times.to_pydatetime(),
                units=attrs['time'].get('units', DEFAULT_TIME_UNITS),
                calendar=time_var.calendar)

            dtype_map = {'float64': 'f8', 'float32': 'f4',
                         'int64': 'i8', 'int32': 'i4'}
            for dim in self.domain['mask'].dims:
                dim_vals = self.domain[dim].values
                dim_dtype = dtype_map.get(
                    str(dim_vals.dtype), self.params['out_precision'])
                dim_var = ncout.createVariable(dim, dim_dtype, (dim, ))
                dim_var[:] = dim_vals

            for p in ['elev', 'lat', 'lon']:
                if p in self.params:
                    self.params.pop(p)
            for k, v in self.params.items():
                # Need to convert some parameters to strings
                if k in ['start', 'stop', 'utc_offset']:
                    v = str(v)
                elif k in ['state_start', 'state_stop', 'out_freq']:
                    # skip
                    continue
                # Don't include complex types
                if isinstance(v, dict):
                    v = json.dumps(v)
                elif not isinstance(v, str) and isinstance(v, Iterable):
                    v = ', '.join(v)

                if isinstance(v, str):
                    v = v.replace("'", "").replace('"', "")
                attrs['_global'][k] = v

            # set global attrs
            for key, val in attrs['_global'].items():
                setattr(ncout, key, val)

            # set variable attrs
            for varname in ncout.variables:
                for key, val in attrs.get(varname, {}).items():
                    setattr(ncout.variables[varname], key, val)

    def write_chunk(self, locks=None):
        '''write data from a single chunk'''
        if not len(self.params['out_vars']):
            return
        for times in self._times:
            filename = self._get_output_filename(times)
            lock = locks.get(filename, DummyLock())
            time_slice = slice(times[0], times[-1])
            with lock:
                with Dataset(filename, mode="r+") as ncout:
                    for varname in self.params['out_vars']:
                        dims = ncout.variables[varname].dimensions[1:]
                        write_slice = ((slice(None), ) + tuple(
                            self._domain_slice[d] for d in dims))
                        ncout.variables[varname][write_slice] = (
                            self.output[varname].sel(time=time_slice).values)

    def run_slice(self):
        """
        Run a single slice of
        """
        self._validate_setup()
        self.disagg = int(self.params['time_step']) < cnst.MIN_PER_DAY
        self.method = MetSim.methods[self.params['method']]
        self.setup_output()
        times = self.met_data['time']
        params = self.params.copy()
        for index, mask_val in np.ndenumerate(self.domain['mask'].values):
            if mask_val > 0:
                locs = {d: i for d, i in zip(self.domain['mask'].dims, index)}
                if self.params['prec_type'].upper() in ['TRIANGLE', 'MIX']:
                    # add variables for triangle precipitation disgregation
                    # method to parameters
                    params['dur'], params['t_pk'] = (
                        add_prec_tri_vars(self.domain.isel(**locs)))
            else:
                continue

            df, state = wrap_run_cell(self.method.run, params,
                                      self.met_data.isel(**locs),
                                      self.state.isel(**locs),
                                      self.disagg, times)

            # Cut the returned data down to the correct time index
            self._unpack_state(state, locs)
            for varname in self.params['out_vars']:
                self.output[varname][locs] = df[varname].values

    def _unpack_state(self, result: pd.DataFrame, locs: dict):
        """Put restart values in the state dataset"""
        # We concatenate with the old state values in case we don't
        # have 90 new days to use
        tmin = np.concatenate((self.state['t_min'].isel(**locs).values[:],
                               result['t_min'].values))
        tmax = np.concatenate((self.state['t_max'].isel(**locs).values[:],
                               result['t_max'].values))
        prec = np.concatenate((self.state['prec'].isel(**locs).values[:],
                               result['prec'].values))
        self.state['t_min'].isel(**locs).values[:] = tmin[-90:]
        self.state['t_max'].isel(**locs).values[:] = tmax[-90:]
        self.state['prec'].isel(**locs).values[:] = prec[-90:]
        state_start = result.index[-1] - pd.Timedelta('89 days')
        self.state['time'].values = date_range(
            state_start, result.index[-1], calendar=self.params['calendar'])

    def _get_output_times(self, freq=None):
        """
        Generate chunked time vectors

        Parameters
        ----------
        freq:
            Output frequency. Given as a Pandas timegrouper string.
            If not given, the entire timeseries will be used.

        Returns
        -------
        times:
            A list of timeseries which represent each of times that
            output files will be created for.
        """
        prototype = self.met_data
        self.disagg = int(self.params['time_step']) < cnst.MIN_PER_DAY

        if self.disagg:
            delta = pd.Timedelta('1 days') - pd.Timedelta(
                '{} minutes'.format(self.params['time_step']))
        else:
            delta = pd.Timedelta('0 days')

        start = pd.Timestamp(prototype['time'].values[0]).to_pydatetime()
        stop = pd.Timestamp(prototype['time'].values[-1]).to_pydatetime()
        times = date_range(start, stop + delta,
                           freq="{}T".format(self.params['time_step']),
                           calendar=self.params['calendar'])

        if freq is None or freq == '':
            times = [times]
        else:
            dummy = pd.Series(np.arange(len(times)), index=times)
            grouper = pd.Grouper(freq=freq)
            times = [t.index for k, t in dummy.groupby(grouper)]
        return times

    def _get_output_filename(self, times):
        suffix = self.get_nc_output_suffix(times)
        fname = '{}_{}.nc'.format(self.params['out_prefix'], suffix)
        output_filename = os.path.join(
            os.path.abspath(self.params['out_dir']), fname)
        return output_filename

    def setup_output(self):

        # output times
        times = self._get_output_times(freq=None)[0]

        # Number of timesteps
        n_ts = len(times)

        shape = (n_ts, ) + self.domain['mask'].shape
        dims = ('time', ) + self.domain['mask'].dims
        coords = {'time': times, **self.domain['mask'].coords}
        self.output = xr.Dataset(coords=coords)
        self.output['time'].encoding['calendar'] = self.params['calendar']

        dtype = self.params['out_precision']
        for varname in self.params['out_vars']:
            self.output[varname] = xr.DataArray(
                data=np.full(shape, np.nan, dtype=dtype),
                coords=coords, dims=dims,
                name=varname, attrs=attrs.get(varname, {}))
        self.output['time'].attrs.update(attrs['time'])

    def _aggregate_state(self):
        """Aggregate data out of the state file and load it into `met_data`"""
        # Precipitation record

        assert self.state.dims['time'] == 90, self.state['time']
        record_dates = date_range(self.params['state_start'],
                                  self.params['state_stop'],
                                  calendar=self.params['calendar'])
        trailing = self.state['prec']
        trailing['time'] = record_dates
        total_precip = xr.concat([trailing, self.met_data['prec']],
                                 dim='time').load()
        total_precip = (cnst.DAYS_PER_YEAR * total_precip.rolling(
            time=90).mean().sel(time=slice(self.params['start'],
                                           self.params['stop'])))

        self.met_data['seasonal_prec'] = total_precip

        # Smoothed daily temperature range
        trailing = self.state['t_max'] - self.state['t_min']

        trailing['time'] = record_dates
        dtr = self.met_data['t_max'] - self.met_data['t_min']
        if (dtr < 0).any():
            raise ValueError("Daily maximum temperature lower"
                             " than daily minimum temperature!")
        sm_dtr = xr.concat([trailing, dtr], dim='time').load()
        sm_dtr = sm_dtr.rolling(time=30).mean().drop(record_dates, dim='time')
        self.met_data['dtr'] = dtr
        self.met_data['smoothed_dtr'] = sm_dtr

    def _validate_setup(self):
        """Updates the global parameters dictionary"""
        errs = [""]

        # Make sure there's some input
        if not len(self.params.get('forcing', [])):
            errs.append("Requires input forcings to be specified")

        # Make sure there is at least one forcing_var
        # They cannot all be constant since we use one as a template
        # for the others
        if not len(self.params.get('forcing_vars', [])):
            errs.append("Requires at least one non-constant forcing")

        # Parameters that can't be empty strings or None
        non_empty = ['out_dir', 'time_step', 'forcing_fmt']
        for each in non_empty:
            if self.params.get(each, None) is None or self.params[each] == '':
                errs.append("Cannot have empty value for {}".format(each))

        # Make sure time step divides evenly into a day
        if (cnst.MIN_PER_DAY % int(self.params.get('time_step', -1)) or
                (int(self.params['time_step']) > (6 * cnst.MIN_PER_HOUR) and
                 int(self.params['time_step']) != cnst.MIN_PER_DAY)):
            errs.append("Time step must be evenly divisible into 1440 "
                        "minutes (24 hours) and less than 360 minutes "
                        "(6 hours). Got {}.".format(self.params['time_step']))

        # Check for required input variable specification
        if self.met_data is not None:
            required_in = ['t_min', 't_max', 'prec']
            for each in required_in:
                if each not in self.met_data.variables:
                    errs.append("Input requires {}".format(each))

        # Make sure that we are going to write out some data
        if not len(self.params.get('out_vars', [])):
            errs.append("Output variable list must not be empty")

        # Check output variables are valid
        daily_out_vars = ['t_min', 't_max', 'prec', 'vapor_pressure',
                          'shortwave', 'tskc', 'pet', 'wind']
        out_var_check = ['temp', 'prec', 'shortwave', 'vapor_pressure',
                         'air_pressure', 'rel_humid', 'spec_humid',
                         'longwave', 'tskc', 'wind']
        if int(self.params.get('time_step', -1)) == 1440:
            out_var_check = daily_out_vars
        for var in self.params.get('out_vars', []):
            if var not in out_var_check:
                errs.append('Cannot output variable {} at timestep {}'.format(
                    var, self.params['time_step']))

        # Check that the parameters specified are available
        opts = {'out_precision': ['f4', 'f8'],
                'lw_cloud': ['default', 'cloud_deardorff'],
                'lw_type': ['default', 'tva', 'anderson',
                            'brutsaert', 'satterlund',
                            'idso', 'prata']}
        for k, v in opts.items():
            if not self.params.get(k, None) in v:
                errs.append("Invalid option given for {}".format(k))

        # If any errors, raise and give a summary
        if len(errs) > 1:
            raise Exception("\n  ".join(errs))

    def get_nc_output_suffix(self, times):
        s, e = times[[0, -1]]
        template = '{:04d}{:02d}{:02d}-{:04d}{:02d}{:02d}'
        return template.format(s.year, s.month, s.day,
                               e.year, e.month, e.day,)


def wrap_run_cell(func: callable, params: dict,
                  ds: xr.Dataset, state: xr.Dataset, disagg: bool,
                  out_times: pd.DatetimeIndex):
    """
    Iterate over a chunk of the domain. This is wrapped
    so we can return a tuple of locs and df.

    Parameters
    ----------
    func: callable
        The function to call to do the work
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

    Returns
    -------
    df_complete: pd.DataFrame
        A dataframe with the disaggregated data in it
    df_base: pd.DataFrame
        A dataframe with the state data in it
    """
    lat = ds['lat'].values.flatten()[0]
    lon = ds['lon'].values.flatten()[0]
    elev = ds['elev'].values.flatten()[0]
    params['elev'] = elev
    params['lat'] = lat
    params['lon'] = lon
    df = ds.to_dataframe()

    # solar_geom returns a tuple due to restrictions of numba
    # for clarity we convert it to a dictionary here
    sg = solar_geom(elev, lat, params['lapse_rate'])
    sg = {'tiny_rad_fract': sg[0], 'daylength': sg[1],
          'potrad': sg[2], 'tt_max0': sg[3]}
    yday = df.index.dayofyear - 1
    df['daylength'] = sg['daylength'][yday]
    df['potrad'] = sg['potrad'][yday]
    df['tt_max'] = sg['tt_max0'][yday]

    # Generate the daily values - these are saved
    # so that we can use a subset of them to write
    # out the state file later
    df_base = func(df, params)

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
            nextday = out_times[-1] + pd.Timedelta('1 days')
            t_end = [ds['t_min'].sel(time=nextday),
                     ds['t_max'].sel(time=nextday)]
        except (KeyError, ValueError):
            # None so that we don't extend the record
            t_end = None

        # Disaggregate to subdaily values
        df_complete = disaggregate(df, params, sg, t_begin, t_end)

        # Calculate the times that we want to get out by chopping
        # off the endpoints that were added on previously
        start = out_times.values[0]
        stop = (out_times.values[-1] + pd.Timedelta('1 days') -
                pd.Timedelta("{} minutes".format(params['time_step'])))
        new_times = date_range(
            start, stop, freq='{}T'.format(params['time_step']),
            calendar=params['calendar'])
    else:
        # convert srad to daily average flux from daytime flux
        df_base['shortwave'] *= df_base['daylength'] / cnst.SEC_PER_DAY
        # If we're outputting daily values, we dont' need to
        # change the output dates - see inside of `if` condition
        # above for more explanation
        new_times = out_times
        df_complete = df_base

    # Cut the returned data down to the correct time index
    df_complete = df_complete.loc[new_times[0]:new_times[-1]]
    return df_complete, df_base


@dask.delayed()
def wrap_run_slice(params, write_locks, domain_slice=NO_SLICE):
    ms = MetSim(params, domain_slice=domain_slice)
    ms.load_inputs()
    ms.run_slice()
    ms.write_chunk(locks=write_locks)


def chunk_domain(chunks, dims):
    '''
    Return a dictionary of chunk slices that can be used to decompose a grid
    '''
    def left(chunk, dim):
        return np.arange(0, dim, chunk)

    def right(chunk, dim):
        nums = np.arange(chunk, dim + chunk, chunk)
        nums[-1] = dim + 1
        return nums

    slices = [[slice(*a) for a in (zip(left(int(chunks[dim]), dims[dim]),
                                       right(int(chunks[dim]), dims[dim])))]
              for dim in chunks.keys()]

    return [dict(zip(chunks.keys(), p)) for p in itertools.product(*slices)]


class DummyLock(object):
    """DummyLock provides the lock API without any actual locking."""
    # This will be available in xarray in the next major version
    def acquire(self, *args):
        pass

    def release(self, *args):
        pass

    def __enter__(self):
        pass

    def __exit__(self, *args):
        pass

    @property
    def locked(self):
        return False


def add_prec_tri_vars(domain):
    """
    Check that variables for triangle precipitation method exist and have
    values that are within allowable ranges. Return these variables.

    Parameters
    ----------
    domain:
        Dataset of domain variables for given location.

    Returns
    -------
    dur
        Array of climatological monthly storm durations. [minutes]
    t_pk
        Array of climatological monthly times to storm peaks. [minutes]
    """
    # Check that variables exist
    try:
        dur = domain['dur']
    except Exception as e:
        raise e("Storm duration and time to peak values are "
                "required in the domain file for the triangle "
                "preciptation disagregation method.")
    try:
        t_pk = domain['t_pk']
    except Exception as e:
        raise e("Storm duration and time to peak values are "
                "required in the domain file for the triangle "
                "preciptation disagregation method.")

    # Check that variable values are within allowable ranges.
    day_length = cnst.MIN_PER_HOUR * cnst.HOURS_PER_DAY
    dur_zero_test = dur <= 0
    dur_day_test = dur > day_length
    if dur_zero_test.any() or dur_day_test.any():
        raise ValueError('Storm duration must be greater than 0 and less than',
                         day_length, '(i.e. the day length in minutes)')

    t_pk_zero_test = t_pk < 0
    t_pk_day_test = t_pk > day_length
    if t_pk_zero_test.any() or t_pk_day_test.any():
        raise ValueError('Storm time to peak must be greater than or equal to '
                         '0, and less than', day_length,
                         '(i.e. the end of a day)')

    return dur, t_pk
