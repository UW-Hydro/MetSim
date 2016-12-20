"""
Disaggregates daily data down to finer grained data using some heuristics
"""

import numpy as np
import pandas as pd
import itertools
import scipy

import metsim
from metsim.configuration import PARAMS as params
from metsim.configuration import CONSTS as consts

def disaggregate(df_daily, solar_geom):
    """
    TODO
    """
    end = metsim.stop + pd.Timedelta('1 days')
    dates_disagg = pd.date_range(metsim.start, end, freq=params['time_step']+'T')
    df_disagg = pd.DataFrame(index=dates_disagg)
    df_disagg['shortwave'] = (shortwave(df_daily['swrad'],
                                       df_daily['dayl'],
                                       df_daily['day_of_year'],
                                       solar_geom['tiny_rad_fract']))
    df_disagg['temp'] = temp(df_daily, df_disagg)
    df_disagg['precip'] = precip(df_daily['precip'])
    df_disagg['longwave'] = longwave(df_daily['lwrad'])
    df_disagg['wind'] = wind(df_daily['wind'])
    return df_disagg

def temp(df_daily, df_disagg):
    """
    TODO
    """
    # Calculate times of min/max temps
    t_Tmin, t_Tmax = set_min_max_hour(df_disagg['shortwave'],
                             len(df_daily['day_of_year']))
    time = np.array(list(next(it) for it in itertools.cycle(
                [iter(t_Tmin), iter(t_Tmax)])))
    temp = np.array(list(next(it) for it in itertools.cycle(
                [iter(df_daily['t_min']), iter(df_daily['t_max'])])))
    cumhours = range(0, 24*len(df_daily['day_of_year']), 24)
    cumhours = np.array(list(next(it) for it in itertools.cycle(
                    [iter(cumhours), iter(cumhours)])))
    time = time + cumhours
    
    interp = scipy.interpolate.PchipInterpolator(time, temp, extrapolate=True)
    temps = interp(range(len(df_disagg.index)))
    return temps


def precip(precip):
    """
    Splits the daily precipitation evenly throughout the day
    """
    return precip.resample(params['time_step']+'T').fillna(method='ffill').sum()


def longwave(longwave):
    """
    Splits the daily longwave evenly throughout the day
    """
    return longwave.resample(params['time_step']+'T').fillna(method='ffill').sum()


def wind(wind):
    """
    Wind is assumed constant throughout the day
    """
    return wind.resample(params['time_step']+'T').fillna(method='ffill').sum()


def shortwave(sw_rad, daylength, day_of_year, tiny_rad_fract):
    """
    TODO
    """
    tiny_step_per_hour = int(consts['SEC_PER_HOUR'] / consts['SRADDT'])
    tmp_rad = sw_rad * daylength / consts['SEC_PER_HOUR'] 
    n_days = len(tmp_rad)
    disaggrad = np.zeros(n_days*consts['HOURS_PER_DAY'] + 1)
    tiny_offset = (params.get("theta_l", 0) - params.get("theta_s", 0) / (consts['HOURS_PER_DAY']/360))
   
    # Tinystep represents a daily set of values - but is constant across days
    tinystep= np.arange(consts['HOURS_PER_DAY'] * tiny_step_per_hour) - tiny_offset
    tinystep[np.array(tinystep<0)] += consts['HOURS_PER_DAY'] * tiny_step_per_hour
    tinystep[np.array(tinystep>(24*tiny_step_per_hour-1))] -= consts['HOURS_PER_DAY'] * tiny_step_per_hour
    chunk_sum = lambda x : np.sum(x.reshape((len(x)/120, 120)), axis=1)
    for day in range(n_days):
        rad = tiny_rad_fract[day_of_year[day]]
        disaggrad[day*consts['HOURS_PER_DAY']: (day+1)*consts['HOURS_PER_DAY']] = (
                chunk_sum(rad[list(tinystep)]) * tmp_rad[day])
    
    # FIXME: Dunno what to do here.
    disaggrad[-2] = disaggrad[-1]
    return disaggrad


def set_min_max_hour(disagg_rad, n_days):
    """
    TODO
    """
    t_Tmax = np.zeros(n_days)
    t_Tmin = np.zeros(n_days)
    for i in range(n_days):
        risehour = sethour = -999
        for hour in range(int(consts['HOURS_PER_DAY']/2)):
            if (disagg_rad[i*consts['HOURS_PER_DAY']+hour] > 0 and
                    (i*consts['HOURS_PER_DAY']+hour==0 or 
                        disagg_rad[int(i*consts['HOURS_PER_DAY'] + hour-1)]<= 0)):
                risehour = hour
        for hour in range(int(consts['HOURS_PER_DAY']/2), int(consts['HOURS_PER_DAY'])):
            if (disagg_rad[i*consts['HOURS_PER_DAY']+hour] <= 0 and 
                    disagg_rad[i*consts['HOURS_PER_DAY']+hour-1]>0):
                sethour = hour
        if i == n_days-1 and sethour == -999:
            sethour = consts['HOURS_PER_DAY'] - 1 
        if risehour >=0 and sethour>=0:
            t_Tmax[i] = 0.67 * (sethour - risehour) + risehour
            t_Tmin[i] = risehour - 1
    return t_Tmin, t_Tmax


