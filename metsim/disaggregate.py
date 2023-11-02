""" Disaggregates daily data down to finer grained data using some heuristics
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
import math
from typing import Tuple

import numpy as np
import pandas as pd
import scipy.interpolate

import metsim.constants as cnst
from metsim.datetime import date_range
from metsim.physics import svp

import sys

def disaggregate(df_daily: pd.DataFrame, params: dict,
                 solar_geom: dict, t_begin: list=None,
                 t_end: list=None) -> pd.DataFrame:
    """
    Take a daily timeseries and scale it down to a finer
    time scale.

    Parameters
    ----------
    df_daily:
        Dataframe containing daily timeseries.
        Should be the result of one of the methods
        provided in the `methods` directory.
    params:
        A dictionary containing the class parameters
        of the MetSim object.
    solar_geom:
        A dictionary of solar geometry variables
    t_begin:
        List of t_min and t_max for day previous to the
        start of `df_daily`. None indicates no extension
        of the record.
    t_end:
        List of t_min and t_max for day after the end
        of `df_daily`. None indicates no extension of
        the record.

    Returns
    -------
    df_disagg:
        A dataframe with sub-daily timeseries.
    """
    # adjust any longitude values to be within [-180, +180] range
    lon_var = params['domain_vars']['lon']
    params[lon_var] = math.remainder(params[lon_var], 360)

    stop = (df_daily.index[-1] + pd.Timedelta('1 days') -
            pd.Timedelta("{} minutes".format(params['time_step'])))
    dates_disagg = date_range(df_daily.index[0], stop,
                              freq='{}T'.format(params['time_step']),
                              calendar=params['calendar'])
    df_disagg = pd.DataFrame(index=dates_disagg)
    n_days = len(df_daily)
    n_disagg = len(df_disagg)
    ts = int(params['time_step'])
    ### assume passthrough implies shortwave was computed using entire day not just daylight like mtclim uses to derive its shortwave
    if (params['method'] == 'passthrough') & (params.get('sw_averaging', '') != 'daylight'):
      df_daily['daylength'] = 86400.0
    df_disagg['shortwave'] = shortwave(df_daily['shortwave'].values,
                                       df_daily['daylength'].values,
                                       df_daily.index.dayofyear,
                                       solar_geom['tiny_rad_fract'],
                                       params)

    t_Tmin, t_Tmax = set_min_max_hour(solar_geom['tiny_rad_fract'],
                                      df_daily.index.dayofyear - 1,
                                      solar_geom['daylength'],
                                      n_days, ts, params)

    df_disagg['temp'] = temp(
        df_daily['t_min'].values, df_daily['t_max'].values,
        n_disagg, t_Tmin, t_Tmax, ts, t_begin, t_end)

    df_disagg['vapor_pressure'] = vapor_pressure(
        df_daily['vapor_pressure'].values, df_disagg['temp'].values,
        t_Tmin, n_disagg, ts)
    df_disagg['vapor_pressure'] = (df_disagg['vapor_pressure']
                                   .fillna(method='ffill')
                                   .fillna(method='bfill'))

    df_disagg['rel_humid'] = relative_humidity(
        df_disagg['vapor_pressure'].values, df_disagg['temp'].values)

    df_disagg['air_pressure'] = pressure(df_disagg['temp'].values,
                                         params['elev'], params['lapse_rate'])

    df_disagg['spec_humid'] = specific_humidity(
        df_disagg['vapor_pressure'].values,
        df_disagg['air_pressure'].values)

    df_disagg['tskc'] = tskc(df_daily['tskc'].values, ts, params)

    if 'longwave' in df_daily:
        daily_lw = df_daily['longwave']
    else:
        daily_lw = None
    df_disagg['longwave'] = longwave(
        df_disagg['temp'].values, df_disagg['vapor_pressure'].values,
        df_disagg['tskc'].values, params, daily_lw)
    df_disagg['prec'] = prec(df_daily['prec'], df_daily['t_min'], ts, params,
                             df_daily.get('t_pk'), df_daily.get('dur'))

    if 'wind' in df_daily:
        df_disagg['wind'] = wind(df_daily['wind'].values, ts, params)

    if params['period_ending']:
        df_disagg.index += pd.Timedelta('{} minutes'.format(params['time_step']))
    return df_disagg.fillna(method='ffill').fillna(method='bfill')


def set_min_max_hour(tiny_rad_fract: np.ndarray, yday: np.ndarray, daylength: np.ndarray,
                     n_days: int, ts: int, params: dict) -> Tuple[np.ndarray]:
    """
    Determine the time at which min and max temp
    is reached for each day.

    Parameters
    ----------
    tiny_rad_fract:
        Array of fraction of shortwave radiation received
        at a shortened timestep. This should be calculated
        by `metsim.physics.solar_geom`.
    yday:
        Array of day of year for each simulated day.
    n_days:
        Number of days in run.
    ts:
        Timestep of run.
    params:
        Dictionary of parameters to use.  Must contain
        `utc_offset` and `tmax_daylength_fraction`.

    Returns
    -------
    (t_t_min, t_t_max):
        A tuple containing 2 timeseries, corresponding
        to time of min and max temp, respectively
    """

    # calculate minute of sunrise and sunset for each day of the year
    rad_mask = 1 * (tiny_rad_fract > 0)
    mask = np.diff(rad_mask)

    # north of the polar circle radiation values of mask </> 0 are eleminated for
    # sunset/sunrise resulting in an array containing less than 365 days for one year
    rise_times = np.zeros(mask.shape[0])
    loc, mult = np.where(mask > 0)
    for i, j in zip(loc, mult):
        rise_times[i] = j * (cnst.SW_RAD_DT / cnst.SEC_PER_MIN)

    set_times = np.zeros(mask.shape[0])
    loc, mult = np.where(mask < 0)
    for i, j in zip(loc, mult):
        set_times[i] = j * (cnst.SW_RAD_DT / cnst.SEC_PER_MIN)

    # Set rise and set times to be equally spaced when in polar region
    day_tot = ((cnst.SEC_PER_DAY / cnst.SW_RAD_DT)
               * (cnst.SW_RAD_DT / cnst.SEC_PER_MIN))
    rise_times[rise_times == 0] = day_tot / 6
    set_times[set_times == 0] = 4 * day_tot / 6

    if params['utc_offset']:
        # rad_fract_per_day = int(cnst.SEC_PER_DAY/cnst.SW_RAD_DT)
        utc_offset = int(((params.get("lon", 0) - params.get("theta_s", 0)) /
                          cnst.DEG_PER_REV) * cnst.MIN_PER_DAY)
        rise_times -= utc_offset
        set_times -= utc_offset

    # map the daily sunrise and sunset to a monotonic timseries (in minutes)
    offset = np.arange(0, n_days * cnst.MIN_PER_HOUR * cnst.HOURS_PER_DAY,
                       cnst.MIN_PER_HOUR * cnst.HOURS_PER_DAY)
    rise_times = rise_times[yday] + offset
    set_times = set_times[yday] + offset

    # time of maximum and minimum temperature calculated thusly
    t_t_max = (params['tmax_daylength_fraction'] * (set_times - rise_times) +
               rise_times) - ts
    t_t_min = rise_times - ts
    return t_t_min, t_t_max


def temp(t_min: np.array, t_max: np.array, out_len: int,
         t_t_min: np.array, t_t_max: np.array, ts: int,
         t_begin: list=None, t_end: list=None) -> np.array:
    """
    Disaggregate temperature using a Hermite polynomial
    interpolation scheme.

    Parameters
    ----------
    t_min:
        Timeseries of daily minimum temperatures.
    t_max:
        Timeseries of daily maximum temperatures.
    out_len:
        Length of the required output vector.
    t_t_min:
        Times at which minimum daily
        temperatures are reached.
    t_t_max:
        Times at which maximum daily
        temperatures are reached.
    ts:
        Timestep for disaggregation
    t_begin:
        List of t_min and t_max for day previous to the
        start of `df_daily`. None indicates no extension
        of the record.
    t_end:
        List of t_min and t_max for day after the end
        of `df_daily`. None indicates no extension of
        the record.

    Returns
    -------
    temps:
        A sub-daily timeseries of temperature.
    """
    # Calculate times of min/max temps by interweaving arrays
    time = np.ravel(np.column_stack((t_t_min, t_t_max)))
    temp = np.ravel(np.column_stack((t_min, t_max)))

    # Account for end points
    ts_ends = cnst.MIN_PER_HOUR * cnst.HOURS_PER_DAY
    time = np.append(np.insert(time, 0, time[0:2] - ts_ends),
                     time[-2:] + ts_ends)

    # If no start or end data is provided to extend the record repeat values
    # This provides space at the ends so that extrapolation doesn't continue
    # in strange ways at the end points
    if t_begin is None:
        t_begin = temp[0:2]
    if t_end is None:
        t_end = temp[-2:]
    temp = np.append(np.insert(temp, 0, t_begin), t_end)

    # Interpolate the values
    interp = scipy.interpolate.PchipInterpolator(time, temp, extrapolate=True)
    temps = interp(ts * np.arange(0, out_len))
    return temps


def prec(prec: pd.Series, t_min: pd.Series, ts: float, params: dict,
         t_pk=None, dur=None):
    """
    Distributes sub-daily precipitation either evenly (uniform) or with a
    triangular (triangle) distribution, depending upon the chosen method.

    Note: The uniform disaggregation returns only through to the beginning of
          the last day. Final values are filled in using a forward fill in the
          top level disaggregate function

    Parameters
    ----------
    prec:
        Daily timeseries of precipitation. [mm]
    t_min:
        Daily timeseries of minimum daily temperature. [C]
    ts:
        Timestep length to disaggregate down to. [minutes]
    params:
        A dictionary of parameters, which contains
        information about which precipitation disaggregation
        method to use.
    month_of_year:
        Timeseries of index of month of year

    Returns
    -------
    prec:
        A sub-daily timeseries of precipitation. [mm]
    """
    def prec_TRIANGLE(daily_prec, t_min, prec_peak, prec_dur, ts, do_mix):
        '''Triangular disaggregation'''
        ts_per_day = int(cnst.HOURS_PER_DAY * cnst.MIN_PER_HOUR)
        n_days = len(daily_prec)
        disagg_prec = np.zeros(int(n_days*ts_per_day))

        # Loop over days
        for i, (t, P) in enumerate(daily_prec.items()):

            if do_mix and t_min[t] < 0:
                prec_day = P * np.ones(ts_per_day)  / ts_per_day
            else:
                times_day = np.arange(ts_per_day)
                prec_day = np.zeros(ts_per_day)

                # Look up climatology
                t_pk = prec_peak[t]
                t_dur = prec_dur[t]

                # Rising and falling time
                t_start = t_pk - 0.5 * t_dur
                t_stop = t_pk + 0.5 * t_dur
                t_plus = times_day[
                        np.logical_and(times_day <= t_pk, times_day >= t_start)]
                t_minus = times_day[
                        np.logical_and(times_day >= t_pk, times_day <= t_stop)]

                # Begin with relative intensity
                prec_day[t_plus] = np.linspace(0, 1.0, len(t_plus))
                prec_day[t_minus] = np.linspace(1.0, 0, len(t_minus))

                # Scale to input precipitation
                prec_day = (P / np.sum(prec_day)) * prec_day
            disagg_prec[i*ts_per_day:(i+1)*ts_per_day] = prec_day
        return disagg_prec

    def prec_UNIFORM(prec: pd.Series, *args):
        return np.repeat(prec.values / cnst.MIN_PER_DAY, cnst.MIN_PER_DAY)

    prec_function = {
        'UNIFORM': prec_UNIFORM,
        'TRIANGLE': prec_TRIANGLE,
        'MIX': prec_TRIANGLE,
    }
    if params['prec_type'].upper() in ['TRIANGLE', 'MIX']:
        if params['prec_type'].upper() == 'MIX':
            do_mix = True
        else:
            do_mix = False
        P_return = prec_TRIANGLE(prec, t_min, t_pk, dur, ts, do_mix)
    else:
        P_return = prec_UNIFORM(prec)

    if params['utc_offset']:
        offset = int(((params["lon"] / cnst.DEG_PER_REV) * cnst.MIN_PER_DAY))
        P_return = np.roll(P_return, -offset)
    P_return = np.sum(P_return.reshape(-1, int(ts)), axis=1).flatten()
    return P_return


def wind(wind: np.array, ts: int, params: dict) -> np.array:
    """
    Wind is assumed constant throughout the day
    Note: this returns only through to the beginning of the
          last day.  Final values are filled in using a
          forward fill in the top level disaggregate function

    Parameters
    ----------
    wind:
        Daily timeseries of wind
    ts:
        Timestep to disaggregate down to

    Returns
    -------
    wind:
        A sub-daily timeseries of wind
    """
    sd_wind = np.repeat(wind, cnst.MIN_PER_DAY)
    if params['utc_offset']:
        offset = int((params["lon"] / cnst.DEG_PER_REV) * cnst.MIN_PER_DAY)
        sd_wind = np.roll(sd_wind, -offset)

    sd_wind = np.mean(sd_wind.reshape(-1, int(ts)), axis=1).flatten()
    return sd_wind


def pressure(temp: np.array, elev: float, lr: float) -> np.array:
    """
    Calculates air pressure.

    Parameters
    ----------
    temp:
        A sub-daily timeseries of temperature
    elev:
        Elevation
    lr:
        Lapse rate

    Returns
    -------
    pressure:
        A sub-daily timeseries of air pressure (kPa)
    """
    temp_corr = cnst.KELVIN + temp + 0.5 * elev * -lr
    ratio = -(elev * cnst.G_STD) / (cnst.R_DRY * temp_corr)
    return cnst.P_STD * np.exp(ratio) / cnst.MBAR_PER_BAR


def specific_humidity(vapor_pressure: np.array,
                      air_pressure: np.array) -> np.array:
    """
    Calculates specific humidity

    Parameters
    ----------
    vapor_pressure:
        A sub-daily timeseries of vapor pressure (Pa)
    air_pressure:
        A sub-daily timeseries of air pressure (Pa)

    Returns
    -------
    spec_humid:
        A sub-daily timeseries of specific humidity
    """
    mix_rat = (cnst.EPS * vapor_pressure) / (air_pressure - vapor_pressure)
    return mix_rat / (1 + mix_rat)


def relative_humidity(vapor_pressure: np.array,
                      temp: np.array) -> np.array:
    """
    Calculate relative humidity from vapor pressure
    and temperature.

    Parameters
    ----------
    vapor_pressure:
        A sub-daily timeseries of vapor pressure
    temp:
        A sub-daily timeseries of temperature

    Returns
    -------
    rh:
        A sub-daily timeseries of relative humidity
    """
    rh = (cnst.MAX_PERCENT * cnst.MBAR_PER_BAR * (vapor_pressure / svp(temp)))
    rh[rh > cnst.MAX_PERCENT] = cnst.MAX_PERCENT
    return rh


def vapor_pressure(vp_daily: np.array, temp: np.array, t_t_min: np.array,
                   n_out: int, ts: int) -> np.array:
    """
    Calculate vapor pressure.  First a linear interpolation
    of the daily values is calculated.  Then this is compared
    to the saturated vapor pressure calculated using the
    disaggregated temperature. When the interpolated vapor
    pressure is greater than the calculated saturated
    vapor pressure, the interpolation is replaced with the
    saturation value.

    Parameters
    ----------
    vp_daily:
        Daily vapor pressure
    temp:
        Sub-daily temperature
    t_t_min:
        Timeseries of minimum daily temperature
    n_out:
        Number of output observations
    ts:
        Timestep to disaggregate down to

    Returns
    -------
    vp:
        A sub-daily timeseries of the vapor pressure
    """
    # Linearly interpolate the values
    interp = scipy.interpolate.interp1d(t_t_min, vp_daily / cnst.MBAR_PER_BAR,
                                        bounds_error=False)
    vp_disagg = interp(ts * np.arange(0, n_out))

    # Account for situations where vapor pressure is higher than
    # saturation point
    vp_sat = svp(temp) / cnst.MBAR_PER_BAR
    vp_disagg = np.where(vp_sat < vp_disagg, vp_sat, vp_disagg)
    return vp_disagg


def longwave(air_temp: np.array, vapor_pressure: np.array,
             tskc: np.array, params: dict, daily_lw=None) -> np.array:
    """
    Calculate longwave. This calculation can be performed
    using a variety of parameterizations for both the
    clear sky and cloud covered emissivity. Options for
    choosing these parameterizations should be passed in
    via the `params` argument.

    For more information about the options provided in this
    function see:

    .. [1] Bohn, T.J., Livneh, B., Oyler, J.W., Running, S.W., Nijssen, B.
      and Lettenmaier, D.P., 2013. Global evaluation of MTCLIM and
      related algorithms for forcing of ecological and hydrological
      models.  Agricultural and forest meteorology, 176, pp.38-49,
      doi:10.1016/j.agrformet.2013.03.003.

    Parameters
    ----------
    air_temp:
        Sub-daily temperature
    vapor_pressure:
        Sub-daily vapor pressure
    tskc:
        Daily cloud fraction
    params:
        A dictionary of parameters, which contains
        information about which emissivity and cloud
        fraction methods to use.

    Returns
    -------
    lwrad:
        A sub-daily timeseries of the longwave radiation
    """
    emissivity_calc = {
        # TVA 1972
        # Tennessee Valley Authority, 1972. Heat and mass transfer between a
        # water surface and the atmosphere. Tennessee Valley Authority, Norris,
        # TN. Laboratory report no. 14. Water resources research report
        # no. 0-6803.
        'DEFAULT': lambda vp: 0.74 + 0.0049 * vp,
        'TVA': lambda vp: 0.74 + 0.0049 * vp,
        # Anderson 1954
        # Anderson, E.R., 1954. Energy budget studies, water loss
        # investigations: lake Hefner studies. U.S. Geol. Surv. Prof. Pap. 269,
        # 71–119 [Available from U.S. Geological Survey, 807 National Center,
        # Reston, VA 20192.].
        'ANDERSON': lambda vp: 0.68 + 0.036 * np.power(vp, 0.5),
        # Brutsaert 1975
        # Brutsaert, W., 1975. On a derivable formula for long-wave radiation
        # from clear skies. Water Resour. Res. 11 (5), 742–744,
        # doi:10.1029/WR011i005p00742.
        'BRUTSAERT': lambda vp: 1.24 * np.power(vp / air_temp, 0.14285714),
        # Satterlund 1979
        # Satterlund, D.R., 1979. An improved equation for estimating long-wave
        # radiation from the atmosphere. Water Resour. Res. 15 (6), 1649–1650,
        # doi:10.1029/WR015i006p01649.
        'SATTERLUND': lambda vp: 1.08 * (
            1 - np.exp(-1 * np.power(vp, (air_temp / 2016)))),
        # Idso 1981
        # Idso, S.B., 1981. A set of equations for full spectrum and 8- to
        # 14-µm and 10.5- to 12.5-µm, thermal radiation from cloudless skies.
        # Water Resour. Res. 17 (2), 295–304, doi:10.1029/WR017i002p00295.
        'IDSO': lambda vp: 0.7 + 5.95e-5 * vp * np.exp(1500 / air_temp),
        # Prata 1996
        # Prata, A.J., 1996. A new long-wave formula for estimating downward
        # clear-sky radiation at the surface. Q. J. R. Meteor. Soc. 122 (533),
        # 1127–1151, doi:10.1002/qj.49712253306.
        'PRATA': lambda vp: (1 - (1 + (46.5 * vp / air_temp)) * np.exp(
            -np.sqrt((1.2 + 3. * (46.5 * vp / air_temp)))))}
    cloud_calc = {
        # TVA 1972 (see above)
        'DEFAULT': lambda emis, tskc: (1.0 + (0.17 * tskc ** 2)) * emis,
        'TVA': lambda emis, tskc: (1.0 + (0.17 * tskc ** 2)) * emis,
        # Deardorff 1978
        # Deardorff, J.W., 1978. Efficient prediction of ground surface
        # temperature and moisture, with an inclusion of a layer of vegetation.
        # J. Geophys. Res. 83 (N64), 1889–1903, doi:10.1029/JC083iC04p01889.
        'CLOUD_DEARDORFF': lambda emis, tskc: tskc + (1 - tskc) * emis
    }
    # Re-index and fill cloud cover, then convert temps to K
    air_temp = air_temp + cnst.KELVIN
    vapor_pressure = vapor_pressure * 10

    # Calculate longwave radiation based on the options
    emiss_func = emissivity_calc[params['lw_type'].upper()]
    emissivity_clear = emiss_func(vapor_pressure)
    emiss_func = cloud_calc[params['lw_cloud'].upper()]
    emissivity = emiss_func(emissivity_clear, tskc)
    lwrad = emissivity * cnst.STEFAN_B * np.power(air_temp, 4)
    if daily_lw is not  None:
        ts = int(params['time_step'])
        ts_per_day = int(cnst.HOURS_PER_DAY * cnst.MIN_PER_HOUR / ts)
        factor = np.mean(lwrad.reshape(-1, ts_per_day), axis=1).flatten()
        factor = daily_lw.values / factor
        factor = np.repeat(factor, ts_per_day)
        lwrad *= factor
    return lwrad


def tskc(tskc: np.array, ts: int, params: dict) -> np.array:
    """
    Disaggregate cloud fraction with uniform interpolation

    Parameters
    ----------
    tskc:
        Daily cloud fraction
    ts:
        Time step to disaggregate to (in minutes)

    Returns
    -------
    tskc:
        Sub-daily timeseries of cloud fraction
    """
    sd_tskc = np.repeat(tskc, cnst.MIN_PER_DAY)
    if params['utc_offset']:
        offset = int((params["lon"] / cnst.DEG_PER_REV) * cnst.MIN_PER_DAY)
        sd_tskc = np.roll(sd_tskc, -offset)

    sd_tskc = np.mean(sd_tskc.reshape(-1, int(ts)), axis=-1).flatten()
    return sd_tskc


def shortwave(sw_rad: np.array, daylength: np.array, day_of_year: np.array,
              tiny_rad_fract: np.array, params: dict) -> np.array:
    """
    Disaggregate shortwave radiation down to a subdaily timeseries.

    Parameters
    ----------
    sw_rad:
        Daily incoming shortwave radiation
    daylength:
        List of daylength time for each day of year
    day_of_year:
        Timeseries of index of days since Jan-1
    tiny_rad_fract:
        Fraction of the daily potential radiation
        during a radiation time step defined by SW_RAD_DT
    params:
        Dictionary of parameters from the MetSim object

    Returns
    -------
    disaggrad:
        A sub-daily timeseries of shortwave radiation.
    """
    ts = int(params['time_step'])
    ts_hourly = float(ts) / cnst.MIN_PER_HOUR
    if params['method'] == 'mtclim' or params.get('sw_averaging', '') == 'daylight':
        tmp_rad = (sw_rad * daylength) / (cnst.SEC_PER_HOUR * ts_hourly)
    else:
        ### if passthrough, still want to do this...but rather than using dailylight daylength uses entire day     
        tmp_rad = (sw_rad * daylength) / (cnst.SEC_PER_HOUR * ts_hourly) #tmp_rad = sw_rad * 24
    n_days = len(tmp_rad)
    ts_per_day = int(cnst.HOURS_PER_DAY * cnst.MIN_PER_HOUR / ts)
    disaggrad = np.zeros(int(n_days * ts_per_day))
    rad_fract_per_day = int(cnst.SEC_PER_DAY / cnst.SW_RAD_DT)
    tmp_rad = np.repeat(tmp_rad, rad_fract_per_day)
    if params['utc_offset']:
        utc_offset = int((params['lon'] / cnst.DEG_PER_REV) * rad_fract_per_day)
        tiny_rad_fract = np.roll(tiny_rad_fract.flatten(), -utc_offset)
        tmp_rad = np.roll(tmp_rad.flatten(), -utc_offset)
        tiny_rad_fract = tiny_rad_fract.flatten()
    else:
        utc_offset = 0
        tiny_rad_fract = tiny_rad_fract.flatten()
    chunk_size = int(ts * (cnst.SEC_PER_MIN / cnst.SW_RAD_DT))
    ts_id = np.repeat(np.arange(ts_per_day), chunk_size)
    for day in range(n_days):
        # Mask to select out from tiny_rad_fract
        radslice = slice((day_of_year[day] - 1) * rad_fract_per_day,
                         (day_of_year[day]) * rad_fract_per_day)
        rad = tiny_rad_fract[radslice]

        # Mask to select out time chunk to place disaggregated values into
        dslice = slice(int(day * ts_per_day), int((day + 1) * ts_per_day))

        # Mask to weight daily solar radiation with
        weight_slice = slice(int(day * rad_fract_per_day),
                             int((day + 1) * rad_fract_per_day))

        rad_chunk = np.bincount(ts_id, weights=rad * tmp_rad[weight_slice])
        disaggrad[dslice] = rad_chunk

    return disaggrad
