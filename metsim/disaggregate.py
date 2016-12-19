"""
Disaggregates daily data down to hourly data using some heuristics
"""

import numpy as np
import pandas as pd

import metsim
from metsim.configuration import PARAMS as params
from metsim.configuration import CONSTS as consts

def disaggregate(df_daily, solar_geom):
    """
    TODO
    """
    end = metsim.stop + pd.Timedelta('1 days')
    dates_hourly = pd.date_range(metsim.start, end, freq=params['time_step']+'T')
    df_disagg = pd.DataFrame(index=dates_hourly)
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
    hours = set_min_max_hour(df_disagg['shortwave'],
                             len(df_daily['day_of_year']))
    df_daily['t_Tmin'] = hours['t_Tmin']
    df_daily['t_Tmax'] = hours['t_Tmax']

    # Fit hermite polynomial and sample daily
    # TODO: Implement this
    # FIXME: This relies on shortwave being implemented correctly
    # x = T_i * day_number
    # y = [Tmin, Tmax, Tmin, Tmax, ... ]
    # interp = scipy.interpolate.PchipInterpolator(x, y, extrapolate=True)
    # temps = interp(range(len(df_disagg.index))
    # return temps

def precip(precip):
    """
    Splits the daily precipitation evenly throughout the day
    """
    return precip.resample(params['time_step']+'T').sum().fillna(method='ffill')


def longwave(longwave):
    """
    Splits the daily longwave evenly throughout the day
    """
    return longwave.resample(params['time_step']+'T').sum().fillna(method='ffill')


def wind(wind):
    """
    Wind is assumed constant throughout the day
    """
    return wind.resample(params['time_step']+'T').sum().fillna(method='ffill')


def shortwave(sw_rad, daylength, day_of_year, tiny_rad_fract):
    """
    TODO
    """
    tiny_step_per_hour = int(consts['SEC_PER_HOUR'] / consts['SRADDT'])
    tmp_rad = sw_rad * daylength / consts['SEC_PER_DAY'] 
    n_days = len(tmp_rad)
    hourlyrad = np.zeros(n_days*consts['HOURS_PER_DAY'] + 1)
    tiny_offset = (params.get("theta_l", 0) - params.get("theta_s", 0) / (consts['HOURS_PER_DAY']/360))
    
    # Tinystep represents a daily set of values - but is constant across days
    tinystep= np.arange(consts['HOURS_PER_DAY'] * tiny_step_per_hour) - tiny_offset
    tinystep[np.array(tinystep<0)] += consts['HOURS_PER_DAY'] * tiny_step_per_hour
    tinystep[np.array(tinystep>(24*tiny_step_per_hour-1))] -= consts['HOURS_PER_DAY'] * tiny_step_per_hour
    chunk_sum = lambda x : np.sum(x.reshape((len(x)/120, 120)), axis=1)
    for day in range(n_days):
        rad = tiny_rad_fract[day_of_year[day]]
        hourlyrad[day*consts['HOURS_PER_DAY']: (day+1)*consts['HOURS_PER_DAY']] = (
                chunk_sum(rad[list(tinystep)]) * tmp_rad[day])
    
    # FIXME: Dunno what to do here.
    hourlyrad[-2] = hourlyrad[-1]
    return hourlyrad


def set_min_max_hour(hourly_rad, n_days):
    """
    TODO
    """
    t_Tmax = np.zeros(n_days)
    t_Tmin = np.zeros(n_days)
    for i in range(n_days):
        risehour = sethour = -999
        for hour in range(int(consts['HOURS_PER_DAY']/2)):
            if (hourly_rad[i*consts['HOURS_PER_DAY']+hour] > 0 and
                    (i*consts['HOURS_PER_DAY']+hour==0 or 
                        hourly_rad[int(i*consts['HOURS_PER_DAY'] + hour-1)]<= 0)):
                risehour = hour
        for hour in range(int(consts['HOURS_PER_DAY']/2), int(consts['HOURS_PER_DAY'])):
            if (hourly_rad[i*consts['HOURS_PER_DAY']+hour] <= 0 and 
                    hourly_rad[i*consts['HOURS_PER_DAY']+hour-1]>0):
                sethour = hour
        if i == n_days -1 and sethour == -999:
            sethour = consts['HOURS_PER_DAY'] - 1 
        if risehour >=0 and sethour>=0:
            t_Tmax[i] = 0.67 * (sethour - risehour) + risehour
            t_Tmin[i] = risehour - 1
    return {'t_Tmin' : t_Tmin, 't_Tmax' : t_Tmax}


