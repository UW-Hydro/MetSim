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

import numpy as np
import pandas as pd
import itertools
import scipy.interpolate

import metsim.constants as cnst
from metsim.physics import svp
from metsim.datetime import date_range


def disaggregate(df_daily: pd.DataFrame, params: dict,
                 solar_geom: dict, t_begin: list=None,
                 t_end: list=None):
    """
    Take a daily timeseries and scale it down to a finer
    time scale.

    Parameters
    ----------
    df_daily: pd.DataFrame
        Dataframe containing daily timeseries.
        Should be the result of one of the methods
        provided in the `methods` directory.
    params: dict
        A dictionary containing the class parameters
        of the MetSim object.
    solar_geom: dict
        A dictionary of solar geometry variables
    t_begin: list
        List of t_min and t_max for day previous to the
        start of `df_daily`. None indicates no extension
        of the record.
    t_end: list
        List of t_min and t_max for day after the end
        of `df_daily`. None indicates no extension of
        the record.

    Returns
    -------
    df_disagg:
        A dataframe with sub-daily timeseries.
    """
    stop = (df_daily.index[-1] + pd.Timedelta('1 days')
            - pd.Timedelta("{} minutes".format(params['time_step'])))
    dates_disagg = date_range(df_daily.index[0], stop,
                              freq='{}T'.format(params['time_step']))
    df_disagg = pd.DataFrame(index=dates_disagg)
    n_days = len(df_daily)
    n_disagg = len(df_disagg)
    ts = float(params['time_step'])
    df_disagg['shortwave'] = shortwave(df_daily['shortwave'],
                                       df_daily['dayl'],
                                       df_daily.index.dayofyear,
                                       solar_geom['tiny_rad_fract'],
                                       params)

    t_Tmin, t_Tmax = set_min_max_hour(df_disagg['shortwave'],
                                      n_days, ts, params)

    df_disagg['temp'] = temp(df_daily, df_disagg, t_Tmin, t_Tmax, ts,
                             t_begin, t_end)

    df_disagg['vapor_pressure'] = vapor_pressure(df_daily['vapor_pressure'],
                                                 df_disagg['temp'],
                                                 t_Tmin, n_disagg, ts)

    df_disagg['rel_humid'] = relative_humidity(df_disagg['vapor_pressure'],
                                               df_disagg['temp'])

    df_disagg['air_pressure'] = pressure(df_disagg['temp'],
                                         params['elev'], params['lapse_rate'])

    df_disagg['spec_humid'] = specific_humidity(df_disagg['vapor_pressure'],
                                                df_disagg['air_pressure'])

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
               rise_times) + ts
    t_t_min = rise_times
    return t_t_min, t_t_max


def temp(df_daily: pd.DataFrame, df_disagg: pd.DataFrame,
         t_t_min: np.array, t_t_max: np.array, ts: float,
         t_begin: list=None, t_end: list=None):
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
    t_begin: list
        List of t_min and t_max for day previous to the
        start of `df_daily`. None indicates no extension
        of the record.
    t_end: list
        List of t_min and t_max for day after the end
        of `df_daily`. None indicates no extension of
        the record.

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
    temps = interp(ts * np.arange(0, len(df_disagg.index)))
    return temps


def prec(prec: pd.Series, ts: float):
    """
    Splits the daily precipitation evenly throughout the day
    Note: this returns only through to the beginning of the
          last day.  Final values are filled in using a
          forward fill in the top level disaggregate function

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
    return wind.resample('{:0.0f}T'.format(ts)).fillna(method='ffill')


def pressure(temp: pd.Series, elev: float, lr: float):
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


def specific_humidity(vapor_pressure: pd.Series, air_pressure: pd.Series):
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
    rh = (cnst.MAX_PERCENT * cnst.MBAR_PER_BAR
          * (vapor_pressure / svp(temp.values)))
    return rh.where(rh < cnst.MAX_PERCENT, cnst.MAX_PERCENT)


def vapor_pressure(vp_daily: pd.Series, temp: pd.Series,
                   t_t_min: np.array, n_out: int, ts: float):
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
    interp = scipy.interpolate.interp1d(t_t_min, vp_daily/cnst.MBAR_PER_BAR,
                                        fill_value='extrapolate')
    vp_disagg = interp(ts * np.arange(0, n_out))

    # Account for situations where vapor pressure is higher than
    # saturation point
    vp_sat = svp(temp.values) / cnst.MBAR_PER_BAR
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
    (lwrad, tskc):
        A sub-daily timeseries of the longwave radiation
        as well as a sub-daily timeseries of the cloud
        cover fraction.
    """
    emissivity_calc = {
        # TVA 1972
        # Tennessee Valley Authority, 1972. Heat and mass transfer between a
        # water surface and the atmosphere. Tennessee Valley Authority, Norris,
        # TN. Laboratory report no. 14. Water resources research report
        # no. 0-6803.
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
        'PRATA': lambda vp: (1 - (1 + (46.5*vp/air_temp)) * np.exp(
            -np.sqrt((1.2 + 3. * (46.5*vp / air_temp)))))
        }
    cloud_calc = {
        # TVA 1972 (see above)
        'TVA': lambda emis: (1.0 + (0.17 * tskc ** 2)) * emis,
        # Deardorff 1978
        # Deardorff, J.W., 1978. Efficient prediction of ground surface
        # temperature and moisture, with an inclusion of a layer of vegetation.
        # J. Geophys. Res. 83 (N64), 1889–1903, doi:10.1029/JC083iC04p01889.
        'CLOUD_DEARDORFF': lambda emis: tskc + (1 - tskc) * emis
        }
    # Re-index and fill cloud cover, then convert temps to K
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
    ts = int(params['time_step'])
    ts_hourly = float(ts) / cnst.MIN_PER_HOUR
    tiny_step_per_hour = cnst.SEC_PER_HOUR / cnst.SW_RAD_DT
    tmp_rad = (sw_rad * daylength) / (cnst.SEC_PER_HOUR * ts_hourly)
    n_days = len(tmp_rad)
    ts_per_day = (cnst.HOURS_PER_DAY * cnst.MIN_PER_HOUR / ts)
    disaggrad = np.zeros(int(n_days*ts_per_day))
    tiny_offset = ((params.get("theta_l", 0) - params.get("theta_s", 0)
                   / (cnst.HOURS_PER_DAY / cnst.DEG_PER_REV)))

    # Tinystep represents a daily set of values - but is constant across days
    tinystep = np.arange(cnst.HOURS_PER_DAY * tiny_step_per_hour) - tiny_offset
    inds = np.asarray(tinystep < 0)
    tinystep[inds] += cnst.HOURS_PER_DAY * tiny_step_per_hour
    inds = np.asarray(tinystep > (cnst.HOURS_PER_DAY * tiny_step_per_hour-1))
    tinystep[inds] -= (cnst.HOURS_PER_DAY * tiny_step_per_hour)
    tinystep = np.asarray(tinystep, dtype=np.int32)

    # Chunk sum takes in the distribution of radiation throughout the day
    # and collapses it into chunks that correspond to the desired timestep
    chunk_size = int(ts * (cnst.SEC_PER_MIN / cnst.SW_RAD_DT))

    # Chunk sum takes in the distribution of radiation throughout the day
    # and collapses it into chunks that correspond to the desired timestep
    def chunk_sum(x):
        return np.sum(x.reshape((int(len(x)/chunk_size), chunk_size)), axis=1)

    for day in range(n_days):
        rad = tiny_rad_fract[day_of_year[day] - 1]
        dslice = slice(int(day * ts_per_day), int((day + 1) * ts_per_day))
        rad_chunk = rad[np.asarray(tinystep, dtype=np.int32)]
        disaggrad[dslice] = chunk_sum(rad_chunk) * tmp_rad[day]
    return disaggrad
