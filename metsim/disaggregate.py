"""
Disaggregates daily data down to finer grained data using some heuristics
"""

import numpy as np
import pandas as pd
import itertools
import scipy

import metsim.constants as cnst
from metsim.physics import svp

def disaggregate(df_daily, params, solar_geom):
    """
    TODO
    """
    stop = params['stop'] + pd.Timedelta('1 days')     
    dates_disagg = pd.date_range(params['start'], stop, freq=params['time_step']+'T')
    df_disagg = pd.DataFrame(index=dates_disagg)
    n_days = len(df_daily['day_of_year'])
    n_disagg = len(df_disagg.index)
    ts = float(params['time_step'])
    
    df_disagg['shortwave'] = (shortwave(df_daily['swrad'],
                                        df_daily['dayl'],
                                        df_daily['day_of_year'],
                                        solar_geom['tiny_rad_fract'],
                                        params))
    
    t_Tmin, t_Tmax = set_min_max_hour(df_disagg['shortwave'], n_days, ts)
    df_disagg['temp'] = temp(df_daily, df_disagg, t_Tmin, t_Tmax, ts)
    
    df_disagg['vapor_pressure'] = vapor_pressure(df_daily['vapor_pressure'],
                                                 df_disagg['temp'],
                                                 t_Tmin, n_disagg, ts)
    
    df_disagg['rel_humid'] = relative_humidity(df_disagg['vapor_pressure'],
                                               df_disagg['temp'])
   
    df_disagg['longwave'], df_disagg['tskc'] = longwave(df_disagg['temp'],
                                                        df_disagg['vapor_pressure'],
                                                        df_daily['tskc'])

    df_disagg['precip'] = precip(df_daily['precip'], ts)
    df_disagg['wind'] = wind(df_daily['wind'], ts)
    df_disagg['pet'] = df_daily['pet'].resample(str(int(ts))+'T').fillna(method='ffill')

    return df_disagg.fillna(method='ffill')


def set_min_max_hour(disagg_rad, n_days, ts):
    """
    TODO
    """
    rad_mask = 1*(disagg_rad > 0)
    diff_mask = np.diff(rad_mask)
    rise_times = np.where(diff_mask>0)[0] * ts 
    set_times = np.where(diff_mask<0)[0] * ts
    t_Tmax = 0.67 * (set_times - rise_times) + rise_times
    t_Tmin = rise_times 
    return t_Tmin, t_Tmax


def temp(df_daily, df_disagg, t_Tmin, t_Tmax, ts):
    """
    TODO
    """
    # Calculate times of min/max temps
    time = np.array(list(next(it) for it in itertools.cycle(
                [iter(t_Tmin), iter(t_Tmax)])))  
    temp = np.array(list(next(it) for it in itertools.cycle(
                [iter(df_daily['t_min']), iter(df_daily['t_max'])])))
    # Account for end points
    ts_ends = cnst.MIN_PER_HOUR * cnst.HOURS_PER_DAY
    time = np.append(np.insert(time, 0, time[0:2]-ts_ends), time[-2:]+ts_ends)
    temp = np.append(np.insert(temp, 0, temp[0:2]), temp[-2:])
    
    # Interpolate the values
    try:
        interp = scipy.interpolate.PchipInterpolator(time, temp, extrapolate=True)
        temps = interp(ts * np.arange(0, len(df_disagg.index)))
    except ValueError:
        temps = np.full(len(df_disagg.index), np.nan)

    return temps


def precip(precip, ts):
    """
    Splits the daily precipitation evenly throughout the day
    """
    scale = int(ts) / (cnst.MIN_PER_HOUR * cnst.HOURS_PER_DAY)
    return (precip*scale).resample(str(int(ts))+'T').fillna(method='ffill')


def wind(wind, ts):
    """
    Wind is assumed constant throughout the day
    """
    return wind.resample(str(int(ts))+'T').fillna(method='ffill')


def relative_humidity(vapor_pressure, temp):
    """
    TODO
    """
    return 100. * np.minimum(vapor_pressure*1000/svp(temp), 1.0)
   

def vapor_pressure(hum_daily, temp, t_Tmin, n_out, ts):
    """Calculate vapor pressure"""
    # Scale down to milibar
    vp_daily = hum_daily / 1000
    # Linearly interpoolate the values
    try:
        interp = scipy.interpolate.interp1d(t_Tmin, vp_daily, fill_value='extrapolate')
        vp_disagg = interp(ts * np.arange(0, n_out))
    except ValueError:
        vp_disagg = np.full(n_out, np.nan)
    # Account for situations where vapor pressure is higher than saturation point
    vp_sat = np.array(svp(temp) / 1000)
    vp_disagg = np.where(vp_sat < vp_disagg, vp_sat, vp_disagg)
    return vp_disagg



def longwave(air_temp, vapor_pressure, tskc):
    """ Calculate longwave """
    emissivity_calc = {
        'DEFAULT'    : lambda x : x,
        'TVA'        : lambda x : 0.74 + 0.0049 * x,
        'ANDERSON'   : lambda x : 0.68 + 0.036 * np.power(x, 0.5),
        'BRUTSAERT'  : lambda x : 1.24 * np.power(x/air_temp, 0.14285714),
        'SATTERLUND' : lambda x : 1.08 * (
            1 - np.exp(-1 * np.power(x, (air_temp/2016)))),
        'IDSO'       : lambda x : 0.7 + 5.95e-5 * x * np.exp(1500/air_temp),
        'PRATA'      : lambda x : (1 - (1 + (46.5*x/air_temp)) *
            np.exp(-np.sqrt((1.2 + 3. * (46.5*x/air_temp)))))
        }
    cloud_calc = {
        'DEFAULT' : lambda x : (1.0 + (0.17 * tskc**2)) * x,
        'CLOUD_DEARDORFF' : lambda x : tskc + (1-tskc)*x
        }
    # Reindex and fill cloud cover, then convert temps to K
    tskc = tskc.reindex(air_temp.index)
    tskc = tskc.fillna(method='ffill')
    air_temp = air_temp + cnst.KELVIN 
    vapor_pressure = vapor_pressure * 10

    # Calculate longwave radiation based on the options 
    emissivity_clear = emissivity_calc[cnst.LW_TYPE.upper()](vapor_pressure)
    emissivity = cloud_calc[cnst.LW_CLOUD.upper()](emissivity_clear) 
    lwrad = emissivity * cnst.STEFAN_B * np.power(air_temp, 4)
    return lwrad, tskc


def shortwave(sw_rad, daylength, day_of_year, tiny_rad_fract, params):
    """
    TODO
    """
    tiny_step_per_hour = cnst.SEC_PER_HOUR / cnst.SRADDT
    tmp_rad = sw_rad * daylength / cnst.SEC_PER_HOUR 
    n_days = len(tmp_rad)
    ts_per_day = cnst.HOURS_PER_DAY * (cnst.MIN_PER_HOUR/int(params['time_step']))
    disaggrad = np.zeros(int(n_days*ts_per_day) + 1)
    tiny_offset = ((params.get("theta_l", 0) - params.get("theta_s", 0) 
                   / (cnst.HOURS_PER_DAY/360)))
    
    # Tinystep represents a daily set of values - but is constant across days
    tinystep= np.arange(cnst.HOURS_PER_DAY * tiny_step_per_hour) - tiny_offset
    tinystep[np.array(tinystep<0)] += cnst.HOURS_PER_DAY * tiny_step_per_hour
    tinystep[np.array(tinystep>(24*tiny_step_per_hour-1))] -= (
            cnst.HOURS_PER_DAY * tiny_step_per_hour)

    # Chunk sum takes in the distribution of radiation throughout the day
    # and collapses it into chunks that correspond to the desired timestep
    chunk_sum = lambda x : np.sum(x.reshape((int(len(x)/120), 120)), axis=1)
    for day in range(n_days):
        rad = tiny_rad_fract[day_of_year[day]]
        disaggrad[int(day*ts_per_day): int((day+1)*ts_per_day)] = (
                chunk_sum(rad[np.array(tinystep).astype(np.int32)]) * tmp_rad[day])
    
    return disaggrad


