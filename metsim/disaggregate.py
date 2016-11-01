"""
Disaggregates daily data down to hourly data using some heuristics
"""

import numpy as np
import pandas as pd

import metsim
from metsim.defaults import PARAMS as params
from metsim.defaults import CONSTS as consts 

tiny_rad_fract = np.zeros(366) #This is updated during the mtclim run

def disaggregate(df_daily):
    """
    TODO
    """
    dates_hourly = pd.date_range(metsim.start, metsim.stop, freq='H') 
    df_hourly = pd.DataFrame(index=dates_hourly)
    _disagg_shortwave(df_daily, df_hourly)
    _disagg_temp(     df_daily, df_hourly)
    _disagg_precip(   df_daily, df_hourly)
    _disagg_thermal(  df_daily, df_hourly)
    _disagg_wind(     df_daily, df_hourly)
    return df_hourly


def _disagg_temp(df_daily, df_hourly):
    """
    TODO
    """
    # Calculate times of min/max temps
    set_min_max_hour(df_daily, df_hourly)
    # Fit hermite polynomial and sample daily 
   


def _disagg_precip(df_daily, df_hourly):
    """
    TODO
    """
    pass


def _disagg_thermal(df_daily, df_hourly):
    """
    TODO
    """
    pass


def _disagg_wind(df_daily, df_hourly):
    """
    TODO
    """   
    pass


def _disagg_shortwave(df_daily, df_hourly):
    """
    TODO
    """
    tiny_step_per_hour = int(3600 / consts['SRADDT'])
    tmp_rad = df_daily['s_swrad']
    n_days = len(tmp_rad)
    hourlyrad = np.zeros(n_days*24+1)
    for i in range(n_days):
        for j in range(24):
            for k in range(tiny_step_per_hour):
                tinystep = j*tiny_step_per_hour + k
                if tinystep < 0:
                    tinystep += 24*tiny_step_per_hour
                if tinystep > 24*tiny_step_per_hour - 1:
                    tinystep -= 24*tiny_step_per_hour
                hourlyrad[i*24+j] += tiny_rad_fract[df_daily['day_of_year'][i]][tinystep]
            #FIXME: This calculation is incorrect
            hourlyrad[i*24+j] *= tmp_rad[i]
    df_hourly['s_swrad'] = hourlyrad


def set_min_max_hour(df_daily, df_hourly):
    """
    TODO
    """   
    hourly_rad = df_hourly['s_swrad']
    n_days = len(df_daily)
    t_max = np.zeros(n_days)
    t_min = np.zeros(n_days)
    for i in range(n_days):
        risehour = sethour = -999
        for hour in range(12):
            if (hourly_rad[i*24+hour] > 0 and 
                    (i*24+hour==0 or hourly_rad[i*24 + hour-1]<= 0)):
                risehour = hour
        for hour in range(12,24):
            if (hourly_rad[i*24+hour] <= 0 and hourly_rad[i*24+hour-1]>0):
                sethour = hour
        if i == n_days -1 and sethour == -999:
            sethour = 23
        if risehour >=0 and sethour>=0:
            t_max[i] - 0.67 * (sethour - risehour) + risehour
            tminhour[i] = rishour - 1
    df_daily['t_Tmin'] = tminhour
    df_daily['t_Tmax'] = tmaxhour

