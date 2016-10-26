"""
Disaggregates daily data down to hourly data using some heuristics
"""

import pandas as pd

import metsim

def disaggregate(df_daily):
    """
    TODO
    """
    dates_hourly = pd.date_range(metsim.start, metsim.stop, freq='H') 
    df_hourly = pd.DataFrame(index=dates_hourly)
    _disagg_temp(   df_daily, df_hourly)
    _disagg_precip( df_daily, df_hourly)
    _disagg_thermal(df_daily, df_hourly)
    _disagg_wind(   df_daily, df_hourly)
    return df_hourly


def _disagg_temp(df_daily, df_hourly):
    """
    TODO
    """
    # Calculate times of min/max temps
    df_daily = pd.concat([df_daily, set_min_max_hour(df_daily, df_hourly)])
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
    n_days = 0
    hourlyrad = np.zeros(n_days*24)
    for i in range(n_days):
        pass
         

def set_min_max_hour(df_daily, df_hourly):
    """
    TODO
    """   
    hourly_rad = df_hourly.get('radiation', None)
    df_daily["t_Tmin"] = 0
    df_daily["t_Tmax"] = 0
    return df_daily


