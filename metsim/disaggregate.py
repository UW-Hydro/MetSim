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
from typing import Tuple

import numpy as np
import pandas as pd
import scipy.interpolate

import metsim.constants as cnst
from metsim.datetime import date_range
from metsim.physics import svp


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
    stop = (df_daily.index[-1] + pd.Timedelta('1 days') -
            pd.Timedelta("{} minutes".format(params['time_step'])))
    dates_disagg = date_range(df_daily.index[0], stop,
                              freq='{}T'.format(params['time_step']))
    df_disagg = pd.DataFrame(index=dates_disagg)
    n_days = len(df_daily)
    n_disagg = len(df_disagg)
    ts = int(params['time_step'])
    df_disagg['shortwave'] = shortwave(df_daily['shortwave'].values,
                                       df_daily['daylength'].values,
                                       df_daily.index.dayofyear,
                                       solar_geom['tiny_rad_fract'],
                                       params)

    t_Tmin, t_Tmax = set_min_max_hour(solar_geom['tiny_rad_fract'],
                                      df_daily.index.dayofyear - 1,
                                      n_days, ts, params)

    df_disagg['temp'] = temp(
        df_daily['t_min'].values, df_daily['t_max'].values,
        n_disagg, t_Tmin, t_Tmax, ts, t_begin, t_end)

    df_disagg['vapor_pressure'] = vapor_pressure(
        df_daily['vapor_pressure'].values, df_disagg['temp'].values,
        t_Tmin, n_disagg, ts)

    df_disagg['rel_humid'] = relative_humidity(
        df_disagg['vapor_pressure'].values, df_disagg['temp'].values)

    df_disagg['air_pressure'] = pressure(df_disagg['temp'].values,
                                         params['elev'], params['lapse_rate'])

    df_disagg['spec_humid'] = specific_humidity(
        df_disagg['vapor_pressure'].values,
        df_disagg['air_pressure'].values)

    df_disagg['tskc'] = tskc(df_daily['tskc'].values, ts, params)

    df_disagg['longwave'] = longwave(
        df_disagg['temp'].values, df_disagg['vapor_pressure'].values,
        df_disagg['tskc'].values, params)

    df_disagg['prec'] = prec(df_daily['prec'].values, df_daily['t_min'].values,
                             ts, params, df_daily.index.month)

    if 'wind' in df_daily:
        df_disagg['wind'] = wind(df_daily['wind'].values, ts, params)

    return df_disagg.fillna(method='ffill')


def set_min_max_hour(tiny_rad_fract: np.array, yday: np.array, n_days: int,
                     ts: int, params: dict) -> Tuple[np.array]:
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
    rise_times = np.where(mask > 0)[1] * (cnst.SW_RAD_DT / cnst.SEC_PER_MIN)
    set_times = np.where(mask < 0)[1] * (cnst.SW_RAD_DT / cnst.SEC_PER_MIN)

    if params['utc_offset']:
        # not used elsewhere:
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
         month_of_year: int):
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
    def prec_TRIANGLE(prec: pd.Series, t_min: pd.Series,
                      month_of_year: int, do_mix: bool, params: dict):

        def P_kernel(t_corners, m, t):
            # calculating precipitation intensity of current time t
            t_start = t_corners[0]
            t_pk = t_corners[1]
            t_end = t_corners[2]
            if (t < t_start) | (t > t_end):
                P = 0
            elif t <= t_pk:
                P = m * (t - t_start)
            elif t >= t_pk:
                P = m * (t_end - t)
            return P

        dur = params['dur']
        t_pk = params['t_pk']
        n_days = len(prec)
        steps_per_day = cnst.MIN_PER_HOUR * cnst.HOURS_PER_DAY
        offset = np.ceil(steps_per_day / 2)
        output_length = int(steps_per_day * n_days)
        index = np.arange(output_length)
        steps_per_two_days = int(((np.ceil(steps_per_day / 2)) * 2) +
                                 steps_per_day)
        # time of step on next day [minutes]
        t_next = 2 * cnst.MIN_PER_DAY + 2
        P_return = pd.Series(np.zeros(output_length, dtype='float'),
                             index=index)

        # create kernel of unit hyetographs, one for each month
        kernels = np.zeros(shape=(12, steps_per_two_days))

        for month in np.arange(12, dtype=int):
            # Calculating the unit hyetograph (kernels) with area of 1 mm

            # Computing monthly constants
            P_pk = 2. / dur[month]  # peak precipitation intensity
            m = P_pk / (dur[month] / 2.)  # slope of precipitation rising limb
            # time of storm start [minutes]
            t_start = (t_pk[month] - (dur[month] / 2.))
            # time of storm end [minutes]
            t_end = (t_pk[month] + (dur[month] / 2.))
            # time step of storm start
            i_start = int(np.floor(t_start) + offset)
            # time step of storm end
            i_end = int(np.floor(t_end) + offset)

            # Initializing timestep variables
            i = i_start  # current timestep
            t = i_start  # current time [minutes]
            t_0 = t     # start time of current timestep [minutes]
            # times of key corners of unit hyetograph
            t_corners = [t_start, t_pk[month], t_end, t_next] + offset
            c = 0   # current corner index
            # end time of current timestep [minutes]
            t_1 = min(t_0, t_corners[c])
            area = 0  # area under curve for current timestep

            # Looping through kernel timesteps
            while i <= i_end:
                P_0 = P_kernel(t_corners, m, t_0)
                P_1 = P_kernel(t_corners, m, t_1)
                area += 0.5 * (P_0 + P_1) * (t_1 - t_0)
                if t_1 == t:  # end of timestep check
                    kernels[month, i] = area
                    area = 0
                    t += 1
                    i += 1
                if t_1 == t_corners[c]:
                    c += 1
                t_0 = t_1
                t_1 = min(t, t_corners[c])

        # Loop through each rain day of the timeseries and apply the kernel for
        # the appropriate month of year
        rain_days = np.asarray(np.where(prec > 0))[0]
        if len(rain_days) > 0:
            for d in rain_days:
                mon = month_of_year[d] - 1
                if do_mix and t_min[d] < 0:
                    i0 = d * steps_per_day
                    i1 = i0 + steps_per_day
                    P_return[i0:i1] += prec[d] * 1 / steps_per_day
                elif d == 0:
                    # beginning of kernel is clipped;
                    # rescale so that clipped kernel sums to original total
                    k0 = int(np.ceil(steps_per_day / 2))
                    i1 = int(np.ceil(steps_per_day * 1.5))
                    P_return[:i1] += (prec[d] / sum(kernels[mon, k0:])
                                      * kernels[mon, k0:])
                elif d == (n_days - 1):
                    # end of kernel is clipped;
                    # rescale so that clipped kernel sums to original total
                    k1 = int(np.ceil(steps_per_day * 1.5))
                    i0 = int(np.floor((d - 0.5) * steps_per_day))
                    P_return[i0:] += (prec[d] / sum(kernels[mon, :k1])
                                      * kernels[mon, :k1])
                else:
                    i0 = int(np.floor((d - 0.5) * steps_per_day))
                    i1 = int(i0 + (2 * steps_per_day))
                    P_return[i0:i1] += prec[d] * kernels[mon]
        P_return = np.around(P_return, decimals=5)
        return P_return.values

    def prec_UNIFORM(prec: pd.Series, *args):
        return np.repeat(prec / cnst.MIN_PER_DAY, cnst.MIN_PER_DAY)

    prec_function = {
        'UNIFORM': prec_UNIFORM,
        'TRIANGLE': prec_TRIANGLE,
        'MIX': prec_TRIANGLE,
    }
    if params['prec_type'].upper() == 'MIX':
        do_mix = True
    else:
        do_mix = False

    P_return = prec_function[params['prec_type'].upper()](prec, t_min,
                                                          month_of_year,
                                                          do_mix, params)

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
        A sub-daily timeseries of vapor pressure (kPa)
    air_pressure:
        A sub-daily timeseries of air pressure (kPa)

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
                                        fill_value='extrapolate')
    vp_disagg = interp(ts * np.arange(0, n_out))

    # Account for situations where vapor pressure is higher than
    # saturation point
    vp_sat = svp(temp) / cnst.MBAR_PER_BAR
    vp_disagg = np.where(vp_sat < vp_disagg, vp_sat, vp_disagg)
    return vp_disagg


def longwave(air_temp: np.array, vapor_pressure: np.array,
             tskc: np.array, params: dict) -> np.array:
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
    tmp_rad = (sw_rad * daylength) / (cnst.SEC_PER_HOUR * ts_hourly)
    n_days = len(tmp_rad)
    ts_per_day = int(cnst.HOURS_PER_DAY * cnst.MIN_PER_HOUR / ts)
    disaggrad = np.zeros(int(n_days * ts_per_day))
    rad_fract_per_day = int(cnst.SEC_PER_DAY / cnst.SW_RAD_DT)
    if params['utc_offset']:
        utc_offset = int(((params.get("lon", 0) - params.get("theta_s", 0)) /
                          cnst.DEG_PER_REV) * rad_fract_per_day)
        tiny_rad_fract = np.roll(tiny_rad_fract.flatten(), -utc_offset)
    else:
        tiny_rad_fract = tiny_rad_fract.flatten()
    chunk_size = int(ts * (cnst.SEC_PER_MIN / cnst.SW_RAD_DT))
    ts_id = np.repeat(np.arange(ts_per_day), chunk_size)
    for day in range(n_days):
        radslice = slice((day_of_year[day] - 1) * rad_fract_per_day,
                         (day_of_year[day]) * rad_fract_per_day)
        rad = tiny_rad_fract[radslice]
        dslice = slice(int(day * ts_per_day), int((day + 1) * ts_per_day))
        rad_chunk = np.bincount(ts_id, weights=rad)
        disaggrad[dslice] = rad_chunk * tmp_rad[day]
    return disaggrad
