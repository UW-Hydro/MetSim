"""
MTCLIM
"""

import numpy as np
import pandas as pd
from warnings import warn

from metsim.configuration import PARAMS  as params
from metsim.configuration import CONSTS  as consts
from metsim.configuration import OPTIONS as options

from metsim.physics import svp, calc_pet, atm_pres

output_variables = ["t_min", "t_max", "precip", "shortwave", "longwave"]

def run(forcing: pd.DataFrame, solar_geom: dict):
    """ 
    Run all of the mtclim forcing generation 
    
    Args:
        forcing: The daily forcings given from input
        solar_geom: Solar geometry of the site
    """
    calc_t_air(forcing)
    calc_precip(forcing)
    calc_snowpack(forcing)
    calc_srad_hum(forcing, solar_geom)
    calc_longwave(forcing)


def calc_t_air(df: pd.DataFrame):
    """ 
    Adjust temperatures according to lapse rates 
    and calculate t_day
    """
    dZ = (params['site_elev'] - params['base_elev'])/1000.
    lapse_rates = [params['t_min_lr'], params['t_max_lr']]
    df['t_max'] = df['t_max'] + dZ * lapse_rates[1]
    df['t_min'] = np.minimum(df['t_min'] + dZ * lapse_rates[0], df['t_max']-0.5)
    t_mean = np.mean(df['t_min'] + df['t_max'])
    df['t_day'] = ((df['t_max'] - t_mean) * params['TDAYCOEF']) + t_mean


def calc_precip(df: pd.DataFrame):
    """ Adjust precipitation according to isoh """
    df['precip'] = (df['precip'] * 
            (df.get('site_isoh', 1) / df.get('base_isoh', 1)))


def calc_snowpack(df: pd.DataFrame):

    def _simple_snowpack(precip, t_min, snowpack=0.0):
        swe = np.array(np.ones(params['n_days']) * snowpack)
        accum = np.array(t_min <= params['SNOW_TCRIT'])
        melt  = np.array(t_min >  params['SNOW_TCRIT'])
        swe[accum] += precip[accum]
        swe[melt]  -= params['SNOW_TRATE'] * (t_min[melt] - params['SNOW_TCRIT'])
        swe = np.maximum(np.cumsum(swe), 0.0) 
        return swe 
   
    df['swe'] = _simple_snowpack(df['precip'], df['t_min'])
    start = (df['day_of_year'] == df['day_of_year'][0])
    end = (df['day_of_year'] == (start-2)%365 + 1) 
    loop_days = np.logical_or(start, end)
    loop_swe  = sum(df['swe'].where(loop_days))/sum(loop_days)

    if np.any(loop_swe):
        df['swe'] = _simple_snowpack(df['precip'], df['t_min'], snowpack=loop_swe)


def calc_srad_hum(df: pd.DataFrame, sg: dict, tol=0.01, win_type='boxcar'):
    """
    Calculate shortwave, humidity  
    """
    def _calc_tfmax(precip, dtr, sm_dtr):
        b = params['B0'] + params['B1'] * np.exp(-params['B2'] * sm_dtr)
        t_fmax = 1.0 - 0.9 * np.exp(-b * np.power(dtr, params['C']))
        inds = np.nonzero(precip > options['SW_PREC_THRESH'])[0]
        t_fmax[inds] *= params['RAIN_SCALAR']
        return t_fmax 

    window = np.zeros(params['n_days'] + 90)
    df['t_max'] = np.maximum(df['t_max'], df['t_min'])
    dtr = df['t_max'] - df['t_min']
    sm_dtr = pd.rolling_window(dtr, window=30, freq='D', 
                               win_type=win_type).fillna(method='bfill')
    if params['n_days'] <= 30:
        warn('Timeseries is shorter than rolling mean window, filling ')
        warn('missing values with unsmoothed data')
        sm_dtr.fillna(dtr, inplace=True)

    sum_precip = df['precip'].values.sum()
    ann_precip = (sum_precip / params['n_days']) * consts['DAYS_PER_YEAR']
    if ann_precip == 0.0:
        ann_precip = 1.0

    if params['n_days'] <= 90:
        sum_precip = df['precip'].values.sum()
        eff_ann_precip = (sum_precip / params['n_days']) * consts['DAYS_PER_YEAR']
        eff_ann_precip = np.maximum(eff_ann_precip, 8.0)
        parray = eff_ann_precip
    else:
        parray = np.zeros(params['n_days'])
        start_yday = df['day_of_year'][0]
        end_yday = df['day_of_year'][-1]
        if start_yday != 1:
            if end_yday == start_yday-1:
                isloop = True
        else:
            if end_yday == 365 or end_yday == 366:
                isloop = True
        
        if isloop:
            for i in range(90):
                window[i] = df['precip'][params['n_days'] - 90 + i]
        else:
            for i in range(90):
                window[i] = df['precip'][i]
        window[90:] = df['precip']

        for i in range(params['n_days']):
            sum_precip = 0.0
            for j in range(90):
                sum_precip += window[i+j]
                sum_precip = (sum_precip / 90.) * consts['DAYS_PER_YEAR']
            sum_precip = np.maximum(sum_precip, 8.0)
            parray[i] = sum_precip

    avg_horizon = (params['site_east_horiz'] + params['site_west_horiz']) / 2.0
    horizon_scalar = 1.0 - np.sin(avg_horizon * consts['RADPERDEG'])
    if (params['site_slope'] > avg_horizon):
        slope_excess = params['site_slope'] - avg_horizon
    else:
        slope_excess = 0.
    if (2.0 * avg_horizon > 180.):
        slope_scalar = 0.
    else:
        slope_scalar = np.clip(1.-(slope_excess / (180.0-2.0 * avg_horizon)), 0, None)
    sky_prop = horizon_scalar * slope_scalar
    df['tfmax'] = _calc_tfmax(df['precip'], dtr, sm_dtr) 

    tdew = df.get('tdew', df['t_min'])
    pva = df.get('hum', svp(tdew))
    pa = atm_pres(params['site_elev'])
    yday = df['day_of_year'] - 1 
    df['dayl'] = sg['daylength'][yday]

    tt_max0 = sg['tt_max0']
    potrad = sg['potrad']
    daylength = sg['daylength']
    yday = df['day_of_year'] - 1
    t_tmax = np.maximum(tt_max0[yday] + (params['ABASE'] * pva), 0.0001)
    t_final = t_tmax * df['tfmax']
    fdir = 1.0 - np.clip(-1.25 * t_final * 1.25, 0., 1.) 
    srad1 = potrad[yday] * t_final * fdir
    srad2 = (potrad[yday] * t_final * (1 - fdir) * \
        (sky_prop + params['DIF_ALB'] * (1.0 - sky_prop)))
    sc = np.zeros_like(df['swe'])
    
    if (options['MTCLIM_SWE_CORR']):
        inds = np.nonzero(df['swe'] > 0. & daylength[yday] > 0.)
        sc[inds] = (1.32 + 0.096 * df['swe'][inds]) *\
            1.0e6 / daylength[yday][inds]
        sc = np.maximum(sc, 100.)  # JJH - this is fishy 
    
    df['swrad'] = srad1 + srad2 + sc
    potrad = (srad1 + srad2 + sc) * daylength[yday] / t_final / consts['SEC_PER_DAY'] 
    tfmax = np.ones(len(sc)) 
    inds = np.nonzero((potrad > 0.) & 
           (df['swrad'] > 0.) & 
           (daylength[yday] > 0))[0]
    tfmax[inds] = np.maximum((df['swrad'][inds] / (potrad[inds] * t_tmax[inds])), 1.)
    df['tfmax'] = tfmax 

    if (options['LW_CLOUD'].upper() == 'CLOUD_DEARDORFF'):
        df['tskc'] = (1. - df['tfmax'])
    else:
        df['tskc'] = np.sqrt((1. - df['tfmax']) / 0.65)

    # Compute PET using SW radiation estimate, and update Tdew, pva **
    tmink = df['t_min'] + consts['KELVIN']
    pet = calc_pet(df['swrad'], df['t_day'], pa, df['dayl'])

    # calculate ratio (PET/effann_prcp) and correct the dewpoint
    ratio = pet / parray
    df['ppratio'] = ratio * 365.25
    tdewk = tmink * (-0.127+1.121 * (1.003-1.444 * ratio+12.312 *
                np.power(ratio, 2)-32.766 * np.power(ratio, 3))+0.0006 * dtr)
    tdew = tdewk - consts['KELVIN']

    pva = svp(tdew)
    if 'hum' not in df:
        df['hum'] = pva
    
    pvs = svp(df['t_day'])
    df['vpd'] = np.maximum(pvs-pva, 0.)


def calc_longwave(df: pd.DataFrame):
    """ Calculate longwave """
    emissivity_calc = {
            'DEFAULT'    : lambda x : x,
            'TVA'        : lambda x : 0.74 + 0.0049 * x,
            'ANDERSON'   : lambda x : 0.68 + 0.036 * np.power(x, 0.5),
            'BRUTSAERT'  : lambda x : 1.24 * np.power(x/air_temp, 0.14285714),
            'SATTERLUND' : lambda x : 1.08 * (1 - np.exp(-1 * np.power(x, (air_temp/2016)))),
            'IDSO'       : lambda x : 0.7 + 5.95e-5 * x * np.exp(1500/air_temp),
            'PRATA'      : lambda x : (1 - (1 + (46.5 * x/air_temp)) *
                np.exp(-1*np.power((1.2 + 3 * (46.5*x/air_temp)), 0.5))),
            }
    cloud_calc = {
            'DEFAULT' : lambda x : (1.0 + (0.17 * tskc**2)) * x,
            'CLOUD_DEARDORFF' : lambda x : tskc * 1.0 + (1-tskc) * x
            }
    tskc = df['tskc']
    air_temp = df['t_day'] + consts['KELVIN'] 
    emissivity_clear = emissivity_calc[options['LW_TYPE'].upper()](air_temp)
    emissivity = cloud_calc[options['LW_CLOUD'].upper()](emissivity_clear) 
    df['lwrad'] = emissivity * consts['STEFAN_B'] * np.power(air_temp, 4)

