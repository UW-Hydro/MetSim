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
    n_disagg = len(df_disagg)
    ts = float(params['time_step'])

    df_disagg['shortwave'] = shortwave(df_daily['swrad'],
                                       df_daily['dayl'],
                                       df_daily['day_of_year'],
                                       solar_geom['tiny_rad_fract'],
                                       params)

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

    return df_disagg.fillna(method='ffill')


def set_min_max_hour(disagg_rad, n_days, ts):
    """
    TODO
    """
    rad_mask = 1*(disagg_rad > 0)
    diff_mask = np.diff(rad_mask)
    rise_times = np.where(diff_mask>0)[0] * ts
    set_times = np.where(diff_mask<0)[0] * ts
    t_t_max = cnst.TMAX_DAYLENGTH_FRACTION * (set_times - rise_times) + rise_times
    t_t_min = rise_times
    return t_t_min, t_t_max


def temp(df_daily, df_disagg, t_t_min, t_t_max, ts):
    """
    TODO
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


def precip(precip, ts):
    """
    Splits the daily precipitation evenly throughout the day
    """
    scale = int(ts) / (cnst.MIN_PER_HOUR * cnst.HOURS_PER_DAY)
    return (precip*scale).resample('{:0.0f}T'.format(ts)).fillna(method='ffill')


def wind(wind, ts):
    """
    Wind is assumed constant throughout the day
    """
    return wind.resample('{:0.0f}T'.format(ts)).fillna(method='ffill')


def relative_humidity(vapor_pressure, temp):
    """
    TODO
    """
    rh = cnst.MAX_PERCENT * cnst.MBAR_PER_BAR * (vapor_pressure/svp(temp))
    return rh.where(rh < cnst.MAX_PERCENT, cnst.MAX_PERCENT) 


def vapor_pressure(vp_daily, temp, t_Tmin, n_out, ts):
    """Calculate vapor pressure"""
    # Linearly interpolate the values
    interp = scipy.interpolate.interp1d(t_Tmin, vp_daily/cnst.MBAR_PER_BAR, 
                                        fill_value='extrapolate')
    vp_disagg = interp(ts * np.arange(0, n_out))

    # Account for situations where vapor pressure is higher than saturation point
    vp_sat = svp(temp) / cnst.MBAR_PER_BAR
    vp_disagg = np.where(vp_sat < vp_disagg, vp_sat, vp_disagg)
    return vp_disagg



def longwave(air_temp, vapor_pressure, tskc):
    """ Calculate longwave """
    emissivity_calc = {
        'DEFAULT'    : lambda vp : vp,
        'TVA'        : lambda vp : 0.74 + 0.0049 * vp,
        'ANDERSON'   : lambda vp : 0.68 + 0.036 * np.power(vp, 0.5),
        'BRUTSAERT'  : lambda vp : 1.24 * np.power(vp/air_temp, 0.14285714),
        'SATTERLUND' : lambda vp : 1.08 * (
            1 - np.exp(-1 * np.power(vp, (air_temp/2016)))),
        'IDSO'       : lambda vp : 0.7 + 5.95e-5 * vp * np.exp(1500/air_temp),
        'PRATA'      : lambda vp : (1 - (1 + (46.5*vp/air_temp)) *
            np.exp(-np.sqrt((1.2 + 3. * (46.5*vp/air_temp)))))
        }
    cloud_calc = {
        'DEFAULT' : lambda emmis : (1.0 + (0.17 * tskc**2)) * emmis,
        'CLOUD_DEARDORFF' : lambda emmis : tskc + (1-tskc)*emmis
        }
    # Reindex and fill cloud cover, then convert temps to K
    tskc = tskc.reindex_like(air_temp).fillna(method='ffill')
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
    tiny_step_per_hour = cnst.SEC_PER_HOUR / cnst.SW_RAD_DT
    tmp_rad = sw_rad * daylength / cnst.SEC_PER_HOUR
    n_days = len(tmp_rad)
    ts_per_day = cnst.HOURS_PER_DAY * (cnst.MIN_PER_HOUR/int(params['time_step']))
    disaggrad = np.zeros(int(n_days*ts_per_day) + 1)
    tiny_offset = ((params.get("theta_l", 0) - params.get("theta_s", 0)
                   / (cnst.HOURS_PER_DAY/cnst.DEG_PER_REV)))

    # Tinystep represents a daily set of values - but is constant across days
    tinystep = np.arange(cnst.HOURS_PER_DAY * tiny_step_per_hour) - tiny_offset
    tinystep[np.array(tinystep<0)] += cnst.HOURS_PER_DAY * tiny_step_per_hour
    tinystep[np.array(tinystep>(cnst.HOURS_PER_DAY*tiny_step_per_hour-1))] -= (
            cnst.HOURS_PER_DAY * tiny_step_per_hour)

    # Chunk sum takes in the distribution of radiation throughout the day
    # and collapses it into chunks that correspond to the desired timestep
    chunk_sum = lambda x : np.sum(x.reshape((int(len(x)/120), 120)), axis=1)
    for day in range(n_days):
        rad = tiny_rad_fract[day_of_year[day]]
        disaggrad[int(day*ts_per_day): int((day+1)*ts_per_day)] = (
                chunk_sum(rad[np.array(tinystep).astype(np.int32)]) * tmp_rad[day])

    return disaggrad


