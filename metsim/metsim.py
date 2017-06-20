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
import struct
import logging
import itertools
import time as tm
from getpass import getuser
from multiprocessing import Pool

import numpy as np
import pandas as pd
import xarray as xr

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

attrs = {'pet': {'units': 'mm d-1', 'long_name': 'potential evaporation',
                 'standard_name': 'water_potential_evaporation_flux'},
         'prec': {'units': 'mm d-1', 'long_name': 'precipitation',
                  'standard_name': 'precipitation_flux'},
         'swrad': {'units': 'W m-2', 'long_name': 'shortwave radiation',
                   'standard_name': 'surface_downwelling_shortwave_flux'},
         'lwrad': {'units': 'W m-2', 'long_name': 'longwave radiation',
                   'standard_name': 'surface_downwelling_longwave_flux'},
         't_max': {'units': 'C', 'long_name': 'maximum daily air temperature',
                   'standard_name': 'daily_maximum_air_temperature'},
         't_min': {'units': 'C', 'long_name': 'minimum daily air temperature',
                   'standard_name': 'daily_minimum_air_temperature'},
         '_global': {'conventions': '1.6', 'title': 'Output from MetSim',
                     'institution': 'University of Washington',
                     'source': 'metsim.py',
                     'history': 'Created: {0} by {1}'.format(now, user),
                     'references': references,
                     'comment': 'no comment at this time'}}


class MetSim(object):
    """
    MetSim handles the distribution of jobs that write to a common file
    by launching muliple processes and queueing up their writeback so that
    work can be done while IO is happening.
    """

    # Class variables
    process_handles = []
    methods = {'mtclim': mtclim}
    params = {
        "method": '',
        "domain": '',
        "state": '',
        "out_dir": '',
        "out_prefix": 'forcing',
        "start": '',
        "stop": '',
        "time_step": '',
        "calendar": 'standard',
        "out_format": '',
        "in_format": None,
        "verbose": 0,
        "base_elev": 0,
        "t_max_lr": '',
        "t_min_lr": '',
        "site_isoh": 1,
        "base_isoh": 1,
        "sw_prec_thresh": 0.0,
        "mtclim_swe_corr": False,
        "lw_cloud": 'cloud_deardorff',
        "lw_type": 'prata',
        "tdew_tol": 1e-6,
        "tmax_daylength_fraction": 0.67,
        "out_vars": ['temp', 'prec', 'shortwave', 'longwave',
                     'vapor_pressure', 'rel_humid']
    }

    def __init__(self, params: dict):
        """
        Constructor
        """
        # Record parameters
        MetSim.params.update(params)
        MetSim.params['dates'] = date_range(params['start'], params['stop'],
                                            calendar=self.params['calendar'])
        logger.setLevel(MetSim.params['verbose'])
        ch.setLevel(MetSim.params['verbose'])
        logger.addHandler(ch)
        self.output = None
        self.met_data = None
        self.ready = False

    def load(self):
        """Load the necessary datasets into memory"""
        # Get the necessary information from the domain
        self.domain = xr.open_dataset(self.params['domain']).rename(
                MetSim.params['domain_vars']).load()
        self.state = xr.open_dataset(self.params['state']).rename(
            MetSim.params.get('state_vars',
                              {'prec': 'prec'})).load()
        self.lat = self.domain['lat']
        self.lon = self.domain['lon']
        self.mask = self.domain['mask']
        self.elev = self.domain['elev']
        self.domain_shape = self.mask.shape
        self.domain_dims = self.mask.dims
        self.disagg = int(self.params['time_step']) < cnst.MIN_PER_DAY
        self.method = MetSim.methods[self.params['method']]
        ilist, jlist = np.nonzero(np.nan_to_num(self.mask.values))

        self.i_idx = ilist
        self.j_idx = jlist

        self.locations = list(zip(ilist, jlist))
        logger.info('found %d locations' % len(self.locations))

        out_dir = MetSim.params['out_dir']
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Input preprocessing
        in_preprocess = {"ascii": self.vic_in_preprocess,
                         "binary": self.vic_in_preprocess,
                         "netcdf": self.netcdf_in_preprocess}
        in_preprocess[MetSim.params['in_format']](self.params['forcing'])
        # Get data from the state file
        self._aggregate_state()
        # Double check that we are ready to do calculations
        self._validate_setup()
        self.ready = True

    def launch(self):
        """Farm out the jobs to separate processes"""
        nprocs = MetSim.params['nprocs']

        # Do the forcing generation and disaggregation if required
        time_dim = pd.to_pydatetime(self.met_data.time.values)
        if self.params['annual']:
            groups = time_dim.groupby(time_dim.year)
        else:
            groups = {'total': time_dim}
        for label, times in groups.items():
            self.pool = Pool(processes=nprocs)
            logger.info("Beginning {}".format(label))
            status = []

            # Add in some end point data for continuity in chunking
            times_ext = times.union([times[0] - pd.Timedelta("1 days"),
                                     times[-1] + pd.Timedelta("1 days")]
                                    ).intersection(time_dim)
            data = self.met_data.sel(time=times_ext)
            self.setup_output(self.met_data.sel(time=times))
            for i, j in self.locations:
                locd = dict(lat=i, lon=j)
                stat = self.pool.apply_async(
                    wrap_run,
                    args=(self.method.run, locd, self.params,
                          data.isel(lat=i, lon=j),
                          self.state.isel(lat=i, lon=j),
                          self.disagg, times, label),
                    callback=self._unpack_results)
                status.append(stat)
            self.pool.close()

            # Check that everything worked
            [stat.get() for stat in status]
            self.pool.join()
            self.write(label)

    def _unpack_results(self, result: tuple):
        """Put results into the master dataset"""
        if len(result) == 3:
            locs, df, state = result
            self._unpack_state(state, locs)
        else:
            locs, df = result
        i = locs['lat']
        j = locs['lon']
        for varname in self.params['out_vars']:
            self.output[varname].values[:, i, j] = df[varname].values

    def _unpack_state(self, result: pd.DataFrame, locs: dict):
        """Put restart values in the state dataset"""
        i = locs['lat']
        j = locs['lon']

        # We concatenate with the old state values in case we don't
        # have 90 new days to use
        tmin = np.concatenate((self.state['t_min'].values[:, i, j],
                               result['t_min'].values))
        tmax = np.concatenate((self.state['t_max'].values[:, i, j],
                               result['t_max'].values))
        prec = np.concatenate((self.state['prec'].values[:, i, j],
                               result['prec'].values))
        self.state['t_min'].values[:, i, j] = tmin[-90:]
        self.state['t_max'].values[:, i, j] = tmax[-90:]
        self.state['prec'].values[:, i, j] = prec[-90:]
        self.state['swe'].values[i, j] = result['swe'].values[-1]
        state_start = result.index[-1] - pd.Timedelta('89 days')
        self.state.time.values = date_range(state_start, result.index[-1],
                                            calendar=self.params['calendar'])

    def run(self):
        """
        Kicks off the disaggregation and queues up data for IO
        """
        time_dim = pd.to_datetime(self.met_data.time.values)
        if self.params['annual']:
            groups = time_dim.groupby(time_dim.year)
        else:
            groups = {'total': time_dim}
        for label, times in groups.items():
            logger.info("Beginning {}".format(label))
            # Add in some end point data for continuity
            times_ext = times.union([times[0] - pd.Timedelta("1 days"),
                                     times[-1] + pd.Timedelta("1 days")]
                                    ).intersection(time_dim)
            data = self.met_data.sel(time=times_ext)
            self.setup_output(self.met_data.sel(time=times))
            for i, j in self.locations:
                locs = dict(lat=i, lon=j)
                logger.info("Processing {}".format(locs))
                ds = data.isel(**locs)
                lat = ds['lat'].values
                elev = ds['elev'].values
                swe = ds['swe'].values
                df = ds.drop(['lat', 'lon', 'elev', 'swe']).to_dataframe()
                # solar_geom returns a tuple due to restrictions of numba
                # for clarity we convert it to a dictionary here
                sg = solar_geom(elev, lat)
                sg = {'tiny_rad_fract': sg[0], 'daylength': sg[1],
                      'potrad': sg[2], 'tt_max0': sg[3]}

                # Generate the daily values - these are saved
                # so that we can use a subset of them to write
                # out the state file later
                df = self.method.run(df, self.params, sg,
                                     elev=elev, swe=swe)

                # Get some values for padding the time list,
                # so that when interpolating in the disaggregation
                # functions we can match endpoints with adjoining
                # chunks - if no data is available, just repeat some
                # default values (this case is used at the very
                # beginning and end of the record)
                if self.disagg:
                    try:
                        prevday = data.time[0] - pd.Timedelta('1 days')
                        t_begin = [self.met_data['t_min'].sel(
                                       time=prevday).isel(lat=i, lon=j),
                                   self.met_data['t_max'].sel(
                                       time=prevday).isel(lat=i, lon=j)]
                    except (KeyError, ValueError):
                        t_begin = [self.state['t_min'].values[-1, i, j],
                                   self.state['t_max'].values[-1, i, j]]
                    try:
                        nextday = pd.datetime(int(label)+1, 1, 1)
                        t_end = [self.met_data['t_min'].sel(
                                     time=nextday).isel(lat=i, lon=j),
                                 self.met_data['t_max'].sel(
                                     time=nextday).isel(lat=i, lon=j)]
                    except (KeyError, ValueError):
                        # None so that we don't extend the record
                        t_end = None

                    self._unpack_state(df, locs)
                    df = disaggregate(df, self.params, sg, t_begin, t_end)
                    start = times[0]
                    stop = (times[-1] + pd.Timedelta('1 days')
                            - pd.Timedelta(self.params['time_step']))
                    new_times = date_range(
                        start, stop,
                        freq='{}T'.format(self.params['time_step']),
                        calendar=self.params['calendar'])
                else:
                    # convert srad to daily average flux from daytime flux
                    self._unpack_state(df, locs)
                    df['swrad'] *= df['dayl'] / cnst.SEC_PER_DAY
                    # If we're outputting daily values, we dont' need to
                    # change the output dates - see inside of `if` condition
                    # above for more explanation
                    new_times = times

                # Cut the returned data down to the correct time index
                self._unpack_results((locs, df.loc[new_times[0]:new_times[-1]]))

            self.write(label)

    def setup_output(self, prototype: xr.Dataset=None):
        if not prototype:
            prototype = self.met_data
        self.output = self.domain.copy(deep=True)
        self.output.attrs = attrs['_global']
        # Number of timesteps
        if self.disagg:
            delta = pd.Timedelta('1 days') - pd.Timedelta(
                '{} minutes'.format(self.params['time_step']))
        else:
            delta = pd.Timedelta('0 days')

        start = pd.Timestamp(prototype.time.values[0]).to_pydatetime()
        stop = pd.Timestamp(prototype.time.values[-1]).to_pydatetime()
        times = date_range(start, stop + delta,
                           freq="{}T".format(MetSim.params['time_step']),
                           calendar=self.params['calendar'])
        n_ts = len(times)

        shape = (n_ts, ) + self.domain_shape
        dims = ('time', ) + self.domain_dims
        coords = {'time': times, **self.domain.mask.coords}
        for varname in MetSim.params['out_vars']:
            self.output[varname] = xr.DataArray(
                data=np.full(shape, np.nan),
                coords=coords, dims=dims,
                name=varname, attrs=attrs.get(varname, {}),
                encoding={'dtype': 'f8', '_FillValue': cnst.FILL_VALUES['f8']})

    def find_elevation(self, lat: float, lon: float) -> float:
        """ Use the domain file to get the elevation """
        return self.elev.sel(lat=lat, lon=lon, method='nearest')

    def _aggregate_state(self):
        """Aggregate data out of the state file and load it into `met_data`"""
        # Precipitation record
        trailing = self.state['prec']
        begin_record = self.params['start'] - pd.Timedelta("90 days")
        end_record = self.params['start'] - pd.Timedelta("1 days")
        record_dates = date_range(begin_record, end_record,
                                  calendar=self.params['calendar'])
        trailing['time'] = record_dates
        total_precip = xr.concat([trailing, self.met_data['prec']], dim='time')
        total_precip = total_precip.rolling(time=90).mean().drop(record_dates,
                                                                 dim='time')
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
        self.met_data['swe'] = xr.Variable(('lat', 'lon'),
                                           self.state['swe'].values)

    def _validate_setup(self):
        """Updates the global parameters dictionary"""
        errs = [""]

        # Make sure there's some input
        if not len(self.params['forcing']):
            errs.append("Requires input forcings to be specified")

        # Parameters that can't be empty strings or None
        non_empty = ['method', 'domain', 'state', 'out_dir',
                     'start', 'stop', 'time_step', 'out_format',
                     'in_format', 't_max_lr', 't_min_lr']
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

        # Convert data types as necessary
        self.params['t_max_lr'] = float(self.params['t_max_lr'])
        self.params['t_min_lr'] = float(self.params['t_min_lr'])

        # Make sure that we are going to write out some data
        if not len(self.params['out_vars']):
            errs.append("Output variable list must not be empty")

        # Check that the parameters specified are available
        opts = {'mtclim_swe_corr': [True, False],
                'lw_cloud': ['default', 'cloud_deardorff'],
                'lw_type': ['default', 'tva', 'anderson',
                            'brutsaert', 'satterlund',
                            'idso', 'prata']}
        for k, v in opts.items():
            if not MetSim.params[k] in v:
                errs.append("Invalid option given for {}".format(k))

        # If any errors, raise and give a summary
        if len(errs) > 1:
            raise Exception("\n  ".join(errs))

    def netcdf_in_preprocess(self, filename: str):
        """Get the extent and load the data"""
        self.met_data = self.read(filename)
        # Subset geographically to match domain
        if len(self.lat) > 0 and len(self.lon) > 0:
            self.met_data = self.met_data.sel(lat=slice(self.lat[0], self.lat[-1]), lon=slice(self.lon[0], self.lon[-1]))

    def vic_in_preprocess(self, job_list: list):
        """Process all files to find spatial extent"""
        # Creates the master dataset which will be used to parallelize
        self.met_data = xr.Dataset(
            coords={'time': MetSim.params['dates'],
                    'lon': self.lon,
                    'lat': self.lat},
            attrs={'n_days': len(MetSim.params['dates'])})
        shape = (len(MetSim.params['dates']), len(self.lat), len(self.lon))

        self.met_data['elev'] = xr.Variable(
            ('lat', 'lon'), np.full((len(self.lat), len(self.lon)), np.nan))
        for var in MetSim.params['in_vars']:
            self.met_data[var] = xr.Variable(
                ('time', 'lat', 'lon'), np.full(shape, np.nan))

        # Fill in the data
        for job in job_list:
            _, lat, lon = os.path.basename(job).split("_")
            i = np.unique(self.i_idx)[list(self.lat.values).index(float(lat))]
            j = np.unique(self.j_idx)[list(self.lon.values).index(float(lon))]
            ds = self.read(job)
            self.met_data['elev'].values[i, j] = self.elev.values[i, j]
            for var in MetSim.params['in_vars']:
                self.met_data[var].values[:, i, j] = ds[var].values

    def write(self, suffix=""):
        """
        Dispatch to the right function based on the configuration given
        """
        dispatch = {
                'netcdf': self.write_netcdf,
                'ascii': self.write_ascii
                }
        self.state.to_netcdf(os.path.join(self.params['out_dir'], 'state.nc'),
                             encoding={'time': {'dtype': 'f8'}})
        dispatch[MetSim.params.get('out_format', 'netcdf').lower()](suffix)

    def write_netcdf(self, suffix: str):
        """Write out as NetCDF to the output file"""
        logger.info("Writing netcdf...")
        fname = '{}_{}.nc'.format(self.params['out_prefix'], suffix)
        output_filename = os.path.join(self.params['out_dir'], fname)
        self.output.to_netcdf(output_filename,
                              unlimited_dims=['time'],
                              encoding={'time': {'dtype': 'f8'}})

    def write_ascii(self, suffix):
        """Write out as ASCII to the output file"""
        logger.info("Writing ascii...")
        shape = self.output.dims
        for i, j in itertools.product(range(shape['lat']),
                                      range(shape['lon'])):
            if self.output.mask[i, j] > 0:
                lat = self.output.lat.values[i]
                lon = self.output.lon.values[j]
                fname = "{}_{}_{}_{}.csv".format(
                            self.params['out_prefix'], lat, lon, suffix)
                fpath = os.path.join(self.params['out_dir'], fname)
                self.output.isel(lat=i, lon=j)[self.params[
                    'out_vars']].to_dataframe().to_csv(fpath)

    def read(self, fpath: str) -> xr.Dataset:
        """
        Dispatch to the right function based on the file extension
        """
        dispatch = {
                'binary': self.read_binary,
                'netcdf': self.read_netcdf,
                'ascii': self.read_ascii
                }
        return dispatch[MetSim.params['in_format']](fpath)

    def read_binary(self, fpath: str) -> xr.Dataset:
        """ Reads a binary forcing file (VIC 4 format) """
        dates = date_range(MetSim.params['start'], MetSim.params['stop'],
                           calendar=self.params['calendar'])
        n_days = len(dates)
        type_lookup = {'signed': 'h', 'unsigned': 'H'}
        # Pack these for nicer syntax in the loop
        var_names = MetSim.params['in_vars'].keys()
        data_list = [[] for var in MetSim.params['in_vars'].keys()]
        n_vars = len(var_names)
        params = [p for p in MetSim.params['in_vars'].values()]
        scales = [float(s.split()[0]) for s in params]
        datatypes = [type_lookup[s.split()[-1]] for s in params]
        with open(fpath, 'rb') as f:
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
        param_list = os.path.basename(fpath).split("_")
        params = {"name": param_list[0],
                  "lat": float(param_list[1]),
                  "lon": float(param_list[2]),
                  "n_days": int(n_days)}
        MetSim.params.update(params)
        params['elev'] = [[self.find_elevation(params['lat'], params['lon'])]]

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

    def read_ascii(self, fpath: str) -> xr.Dataset:
        """Read in an ascii forcing file"""
        dates = date_range(MetSim.params['start'], MetSim.params['stop'],
                           calendar=self.params['calendar'])
        names = MetSim.params['in_vars'].keys()
        ds = pd.read_table(fpath, header=None, delim_whitespace=True,
                           names=names).head(len(dates))
        ds.index = dates
        return ds

    def read_netcdf(self, fpath: str) -> xr.Dataset:
        """
        Read in a NetCDF file and add elevation information
        """
        ds = xr.open_dataset(fpath)
        ds.rename(MetSim.params['in_vars'], inplace=True)
        varlist = list(MetSim.params['in_vars'].values())
        ds = ds[varlist].sel(
            time=slice(MetSim.params['start'], MetSim.params['stop']))
        dates = ds.indexes['time']
        # Add elevation and day of year data
        ds['elev'] = self.elev
        ds['day_of_year'] = xr.Variable(('time', ), dates.dayofyear)

        # Update the configuration
        MetSim.params.update({"n_days": len(ds['time'])})
        return ds.load()


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
    df = ds.drop(['lat', 'lon', 'elev', 'swe']).to_dataframe()
    # solar_geom returns a tuple due to restrictions of numba
    # for clarity we convert it to a dictionary here
    sg = solar_geom(elev, lat)
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
        stop = out_times[-1] + pd.Timedelta('23 hours')
        new_times = date_range(
            start, stop, freq='{}T'.format(params['time_step']),
            calendar=params['calendar'])
    else:
        # convert srad to daily average flux from daytime flux
        df_base['swrad'] *= df_base['dayl'] / cnst.SEC_PER_DAY
        # If we're outputting daily values, we dont' need to
        # change the output dates - see inside of `if` condition
        # above for more explanation
        new_times = out_times

    # Cut the returned data down to the correct time index
    df_complete = df_complete.loc[new_times[0]:new_times[-1]]
    df_base = df_base.loc[new_times[0]:new_times[-1]]
    return (loc, df_complete, df_base)
