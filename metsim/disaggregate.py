"""
Disaggregates daily data down to hourly data using some heuristics
"""

import numpy as np
import pandas as pd

import metsim
from metsim.configuration import PARAMS as params
from metsim.configuration import CONSTS as consts 

#FIXME: These should be initialized before forcings are generated
tiny_rad_fract = np.zeros(366) #This is updated during the mtclim run
daylength      = np.zeros(366) #This is updated during the mtclim run
slope_potrad   = np.zeros(366) #This is updated during the mtclim run
flat_potrad    = np.zeros(366) #This is updated during the mtclim run
tt_max0        = np.zeros(366) #This is updated during the mtclim run

def disaggregate(df_daily):
    """
    TODO
    """
    dates_hourly = pd.date_range(metsim.start, metsim.stop, freq='H') 
    df_hourly = pd.DataFrame(index=dates_hourly)
    print(len(dates_hourly))
    df_hourly['shortwave'] = shortwave(df_daily['s_swrad'], df_daily['day_of_year'])
    df_hourly['temp'] = temp(df_daily, df_hourly)
    df_hourly['precip'] = precip(df_daily['precip'])
    df_hourly['thermal'] = thermal(df_daily, df_hourly)
    df_hourly['wind'] = wind(df_daily['wind'])
    return df_hourly


def temp(df_daily, df_hourly):
    """
    TODO
    """
    # Calculate times of min/max temps
    hours = set_min_max_hour(df_hourly['shortwave'], 
                             len(df_hourly['day_of_year']))
    df_daily['t_Tmin'] = hours['t_Tmin']
    df_daily['t_Tmax'] = hours['t_Tmax']

    # Fit hermite polynomial and sample daily 
    # TODO 


def precip(precip):
    """
    Splits the daily precipitation evenly throughout the day 
    """
    return precip.resample('H', how='sum').fillna(method='ffill')/24. 


def thermal(df_daily, df_hourly):
    """
    TODO
    """
    pass


def wind(wind):
    """
    Wind is assumed constant throughout the day
    """   
    return wind.resample('H', how='sum').fillna(method='ffill')


def shortwave(sw_rad, day_of_year):
    """
    TODO
    """
    tiny_step_per_hour = int(3600 / consts['SRADDT'])
    tmp_rad = sw_rad 
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
                hourlyrad[i*24+j] += tiny_rad_fract[day_of_year[i]][tinystep]
            #FIXME: This calculation is incorrect
            hourlyrad[i*24+j] *= tmp_rad[i]
    print(len(hourlyrad))
    return hourlyrad


def set_min_max_hour(hourly_rad, n_days):
    """
    TODO
    """   
    t_Tmax = np.zeros(n_days)
    t_Tmin = np.zeros(n_days)
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
            t_Tmax[i] - 0.67 * (sethour - risehour) + risehour
            t_Tmin[i] = rishour - 1
    return {'t_Tmin' : t_Tmin, 't_Tmax' : t_Tmax} 


