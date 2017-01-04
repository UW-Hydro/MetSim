"""
Disaggregates daily data down to finer grained data using some heuristics
"""

import numpy as np
import pandas as pd
import itertools
import scipy

from metsim.configuration import PARAMS as params
from metsim.configuration import CONSTS as consts

def disaggregate(df_daily, solar_geom):
    """
    TODO
    """
    stop = params['stop'] + pd.Timedelta('1 days')     
    dates_disagg = pd.date_range(params['start'], stop, freq=params['time_step']+'T')
    df_disagg = pd.DataFrame(index=dates_disagg)
    df_disagg['shortwave'] = (shortwave(df_daily['swrad'],
                                       df_daily['dayl'],
                                       df_daily['day_of_year'],
                                       solar_geom['tiny_rad_fract']))
    df_disagg['temp'] = temp(df_daily, df_disagg)
    df_disagg['precip'] = precip(df_daily['precip'])
    df_disagg['longwave'] = longwave(df_daily['lwrad'])
    df_disagg['wind'] = wind(df_daily['wind'])
    return df_disagg.fillna(method='ffill')


def set_min_max_hour(disagg_rad, n_days):
    """
    TODO
    """
    ts = float(params['time_step'])
    rad_mask = 1*(disagg_rad > 0)
    diff_mask = np.diff(rad_mask)
    rise_times = np.where(diff_mask>0)[0] * ts 
    set_times = np.where(diff_mask<0)[0] * ts
    t_Tmax = 0.67 * (set_times - rise_times) + rise_times
    t_Tmin = rise_times - float(params['time_step'])
    return t_Tmin, t_Tmax


def temp(df_daily, df_disagg):
    """
    TODO
    """
    # Calculate times of min/max temps
    n_days = len(df_daily['day_of_year'])
    t_Tmin, t_Tmax = set_min_max_hour(df_disagg['shortwave'], n_days)
    time = np.array(list(next(it) for it in itertools.cycle(
                [iter(t_Tmin), iter(t_Tmax)])))
    temp = np.array(list(next(it) for it in itertools.cycle(
                [iter(df_daily['t_min']), iter(df_daily['t_max'])])))
    
    try:
        interp = scipy.interpolate.PchipInterpolator(time, temp, extrapolate=True)
        temps = interp(float(params['time_step']) * np.arange(0, len(df_disagg.index)))
    except ValueError:
        temps = np.full(len(df_disagg.index), np.nan)

    return temps


def precip(precip):
    """
    Splits the daily precipitation evenly throughout the day
    """
    scale = int(params['time_step']) / (consts['MIN_PER_HOUR'] * consts['HOURS_PER_DAY'])
    return (precip*scale).resample(params['time_step']+'T').fillna(method='ffill')


def longwave(longwave):
    """
    Splits the daily longwave evenly throughout the day
    """
    return longwave.resample(params['time_step']+'T').fillna(method='ffill')


def wind(wind):
    """
    Wind is assumed constant throughout the day
    """
    return wind.resample(params['time_step']+'T').fillna(method='ffill')


def shortwave(sw_rad, daylength, day_of_year, tiny_rad_fract):
    """
    TODO
    """
    tiny_step_per_hour = int(consts['SEC_PER_HOUR'] / consts['SRADDT'])
    tmp_rad = sw_rad * daylength / consts['SEC_PER_HOUR'] 
    n_days = len(tmp_rad)
    ts_per_day = consts['HOURS_PER_DAY'] * (consts['MIN_PER_HOUR']/int(params['time_step']))
    disaggrad = np.zeros(n_days*ts_per_day + 1)
    tiny_offset = (params.get("theta_l", 0) - params.get("theta_s", 0) / (consts['HOURS_PER_DAY']/360))
   
    # Tinystep represents a daily set of values - but is constant across days
    tinystep= np.arange(consts['HOURS_PER_DAY'] * tiny_step_per_hour) - tiny_offset
    tinystep[np.array(tinystep<0)] += consts['HOURS_PER_DAY'] * tiny_step_per_hour
    tinystep[np.array(tinystep>(24*tiny_step_per_hour-1))] -= consts['HOURS_PER_DAY'] * tiny_step_per_hour
    chunk_sum = lambda x : np.sum(x.reshape((len(x)/120, 120)), axis=1)
    for day in range(n_days):
        rad = tiny_rad_fract[day_of_year[day]]
        disaggrad[day*ts_per_day: (day+1)*ts_per_day] = (
                chunk_sum(rad[list(tinystep)]) * tmp_rad[day])
    
    # FIXME: Dunno what to do here.
    disaggrad[-2] = disaggrad[-1]
    return disaggrad


