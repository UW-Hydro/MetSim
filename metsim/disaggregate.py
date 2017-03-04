"""
Disaggregates daily data down to finer grained data using some heuristics
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

import numpy as np
import pandas as pd
import itertools
import scipy

import metsim.constants as cnst
from metsim.physics import svp


def disaggregate(df_daily: pd.DataFrame, params: dict,
                 solar_geom: dict):
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

    Returns
    -------
    df_disagg:
        A dataframe with sub-daily timeseries.
    """
    stop = params['stop'] + pd.Timedelta('1 days')
    dates_disagg = pd.date_range(params['start'], stop,
                                 freq='{}T'.format(params['time_step']))
    df_disagg = pd.DataFrame(index=dates_disagg)
    n_days = len(df_daily)
    n_disagg = len(df_disagg)
    ts = float(params['time_step'])

    df_disagg['shortwave'] = shortwave(df_daily['swrad'],
                                       df_daily['dayl'],
                                       df_daily.index.dayofyear,
                                       solar_geom['tiny_rad_fract'],
                                       params)

    t_Tmin, t_Tmax = set_min_max_hour(df_disagg['shortwave'],
                                      n_days, ts, params)

    df_disagg['temp'] = temp(df_daily, df_disagg, t_Tmin, t_Tmax, ts)

    df_disagg['vapor_pressure'] = vapor_pressure(df_daily['vapor_pressure'],
                                                 df_disagg['temp'],
                                                 t_Tmin, n_disagg, ts)

    df_disagg['rel_humid'] = relative_humidity(df_disagg['vapor_pressure'],
                                               df_disagg['temp'])

    df_disagg['longwave'], df_disagg['tskc'] = longwave(
        df_disagg['temp'], df_disagg['vapor_pressure'],
        df_daily['tskc'], params)

    df_disagg['prec'] = prec(df_daily['prec'], ts)
    if 'wind' in df_daily:
        df_disagg['wind'] = wind(df_daily['wind'], ts)

    return df_disagg.fillna(method='ffill')


def set_min_max_hour(disagg_rad: pd.Series, n_days: int,
                     ts: float, params: dict):
    """
    Determine the time at which min and max temp
    is reached for each day.

    Parameters
    ----------
    disagg_rad:
        Shortwave radiation disaggregated
        to sub-daily timesteps.
    n_days:
        The number of days being disaggregated
    ts:
        Timestep used for disaggregation
    params:
        A dictionary of class parameters of
        the MetSim object.

    Returns
    -------
    (t_t_min, t_t_max):
        A tuple containing 2 timeseries, corresponding
        to time of min and max temp, respectively
    """
    rad_mask = 1*(disagg_rad > 0)
    diff_mask = np.diff(rad_mask)
    rise_times = np.where(diff_mask > 0)[0] * ts
    set_times = np.where(diff_mask < 0)[0] * ts
    t_t_max = (params['tmax_daylength_fraction'] * (set_times - rise_times) +
               rise_times)
    t_t_min = rise_times
    return t_t_min, t_t_max


def temp(df_daily: pd.DataFrame, df_disagg: pd.DataFrame,
         t_t_min: np.array, t_t_max: np.array, ts: float):
    """
    Disaggregate temperature using a Hermite polynomial
    interpolation scheme.

    Parameters
    ----------
    df_daily:
        A dataframe of daily values.
    df_disagg:
        A dataframe of sub-daily values.
    t_t_min:
        Times at which minimum daily
        temperatures are reached.
    t_t_max:
        Times at which maximum daily
        temperatures are reached.
    ts:
        Timestep for disaggregation

    Returns
    -------
    temps:
        A sub-daily timeseries of temperature.
    """
    # Calculate times of min/max temps
    time = np.array(list(next(it) for it in itertools.cycle(
                [iter(t_t_min), iter(t_t_max)])))
    temp = np.array(list(next(it) for it in itertools.cycle(
                [iter(df_daily['t_min']), iter(df_daily['t_max'])])))
    # Account for end points
    ts_ends = cnst.MIN_PER_HOUR * cnst.HOURS_PER_DAY
    time = np.append(np.insert(time, 0, time[0:2]-ts_ends), time[-2:]+ts_ends)
    temp = np.append(np.insert(temp, 0, temp[0:2]), temp[-2:])

    # Interpolate the values
    interp = scipy.interpolate.PchipInterpolator(time, temp, extrapolate=True)
    temps = interp(ts * np.arange(0, len(df_disagg.index)))
    return temps


def prec(prec: pd.Series, ts: float):
    """
    Splits the daily precipitation evenly throughout the day

    Parameters
    ----------
    prec:
        Daily timeseries of precipitation
    ts:
        Timestep to disaggregate down to

    Returns
    -------
    prec:
        A sub-daily timeseries of precipitation
    """
    scale = int(ts) / (cnst.MIN_PER_HOUR * cnst.HOURS_PER_DAY)
    return (prec * scale).resample(
        '{:0.0f}T'.format(ts)).fillna(method='ffill')


def wind(wind: pd.Series, ts: float):
    """
    Wind is assumed constant throughout the day

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
    return wind.resample('{:0.0f}T'.format(ts)).fillna(method='ffill')


def relative_humidity(vapor_pressure: pd.Series, temp: pd.Series):
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
    rh = cnst.MAX_PERCENT * cnst.MBAR_PER_BAR * (vapor_pressure / svp(temp))
    return rh.where(rh < cnst.MAX_PERCENT, cnst.MAX_PERCENT)


def vapor_pressure(vp_daily: pd.Series, temp: pd.Series,
                   t_t_min: np.array, n_out: int, ts: float):
    """
    Calculate vapor pressure.  First a linear inerpolation
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
    interp = scipy.interpolate.interp1d(t_t_min, vp_daily/cnst.MBAR_PER_BAR,
                                        fill_value='extrapolate')
    vp_disagg = interp(ts * np.arange(0, n_out))

    # Account for situations where vapor pressure is higher than
    # saturation point
    vp_sat = svp(temp) / cnst.MBAR_PER_BAR
    vp_disagg = np.where(vp_sat < vp_disagg, vp_sat, vp_disagg)
    return vp_disagg


def longwave(air_temp: pd.Series, vapor_pressure: pd.Series,
             tskc: pd.Series, params: dict):
    """
    Calculate longwave. This calculation can be performed
    using a variety of parameterizations for both the
    clear sky and cloud covered emissivity. Options for
    choosing these parameterizations should be passed in
    via the `params` argument.

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
    (lwrad, tskc):
        A sub-daily timeseries of the longwave radiation
        as well as a sub-daily timeseries of the cloud
        cover fraction.
    """
    emissivity_calc = {
        'DEFAULT': lambda vp: vp,
        'TVA': lambda vp: 0.74 + 0.0049 * vp,
        'ANDERSON': lambda vp: 0.68 + 0.036 * np.power(vp, 0.5),
        'BRUTSAERT': lambda vp: 1.24 * np.power(vp / air_temp, 0.14285714),
        'SATTERLUND': lambda vp: 1.08 * (
            1 - np.exp(-1 * np.power(vp, (air_temp / 2016)))),
        'IDSO': lambda vp: 0.7 + 5.95e-5 * vp * np.exp(1500 / air_temp),
        'PRATA': lambda vp: (1 - (1 + (46.5*vp/air_temp)) * np.exp(
            -np.sqrt((1.2 + 3. * (46.5*vp / air_temp)))))
        }
    cloud_calc = {
        'DEFAULT': lambda emis: (1.0 + (0.17 * tskc ** 2)) * emis,
        'CLOUD_DEARDORFF': lambda emis: tskc + (1 - tskc) * emis
        }
    # Reindex and fill cloud cover, then convert temps to K
    tskc = tskc.reindex_like(air_temp).fillna(method='ffill')
    air_temp = air_temp + cnst.KELVIN
    vapor_pressure = vapor_pressure * 10

    # Calculate longwave radiation based on the options
    emiss_func = emissivity_calc[params['lw_type'].upper()]
    emissivity_clear = emiss_func(vapor_pressure)
    emiss_func = cloud_calc[params['lw_cloud'].upper()]
    emissivity = emiss_func(emissivity_clear)
    lwrad = emissivity * cnst.STEFAN_B * np.power(air_temp, 4)
    return lwrad, tskc


def shortwave(sw_rad: pd.Series, daylength: pd.Series, day_of_year: pd.Series,
              tiny_rad_fract: np.array, params: dict):
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
    tiny_step_per_hour = cnst.SEC_PER_HOUR / cnst.SW_RAD_DT
    tmp_rad = sw_rad * daylength / cnst.SEC_PER_HOUR
    n_days = len(tmp_rad)
    ts_per_day = (cnst.HOURS_PER_DAY *
                  cnst.MIN_PER_HOUR / int(params['time_step']))
    disaggrad = np.zeros(int(n_days*ts_per_day) + 1)
    tiny_offset = ((params.get("theta_l", 0) - params.get("theta_s", 0) /
                    (cnst.HOURS_PER_DAY / cnst.DEG_PER_REV)))

    # Tinystep represents a daily set of values - but is constant across days
    tinystep = np.arange(cnst.HOURS_PER_DAY * tiny_step_per_hour) - tiny_offset
    inds = np.array(tinystep < 0)
    tinystep[inds] += cnst.HOURS_PER_DAY * tiny_step_per_hour
    inds = np.array(tinystep > (cnst.HOURS_PER_DAY * tiny_step_per_hour-1))
    tinystep[inds] -= (cnst.HOURS_PER_DAY * tiny_step_per_hour)

    # Chunk sum takes in the distribution of radiation throughout the day
    # and collapses it into chunks that correspond to the desired timestep
    def chunk_sum(x):
        return np.sum(x.reshape((int(len(x)/120), 120)), axis=1)

    for day in range(n_days):
        rad = tiny_rad_fract[day_of_year[day] - 1]
        dslice = slice(int(day * ts_per_day), int((day + 1) * ts_per_day))
        disaggrad[dslice] = (
            chunk_sum(rad[np.array(tinystep).astype(np.int32)]) * tmp_rad[day])

    return disaggrad
