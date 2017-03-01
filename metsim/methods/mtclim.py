"""
MTCLIM
"""

import numpy as np
import pandas as pd
from warnings import warn

import metsim.constants as cnst
from metsim.disaggregate import disaggregate
from metsim.physics import svp, calc_pet, atm_pres, solar_geom


def run(forcing: pd.DataFrame, params: dict, elev: float, lat: float,
        disagg=True):
    """
    Run all of the mtclim forcing generation

    Args:
        forcing: The daily forcings given from input
        solar_geom: Solar geometry of the site
    """

    # solar_geom returns a tuple due to restrictions of numba
    # for clarity we convert it to a dictionary here
    sg = solar_geom(elev, lat)
    sg = {'tiny_rad_fract': sg[0],
          'daylength': sg[1],
          'potrad': sg[2],
          'tt_max0': sg[3]}
    params['n_days'] = len(forcing)
    calc_t_air(forcing, elev, params)
    calc_prec(forcing, params)
    calc_snowpack(forcing, params)
    calc_srad_hum(forcing, sg, elev, params)

    if disagg:
        forcing = disaggregate(forcing, params, sg)
    else:
        # convert srad to daily average flux from daytime flux
        print(forcing.head())
        forcing['swrad'] *= forcing['dayl'] / cnst.SEC_PER_DAY

    return forcing


def calc_t_air(df: pd.DataFrame, elev: float, params: dict):
    """
    Adjust temperatures according to lapse rates
    and calculate t_day
    """
    dZ = (elev - params['base_elev'])/cnst.M_PER_KM
    lapse_rates = [params['t_min_lr'], params['t_max_lr']]
    t_max = df['t_max'] + dZ * lapse_rates[1]
    t_min = df['t_min'].where(df['t_min'] + dZ * lapse_rates[0] < t_max-0.5,
                              t_max - 0.5)
    t_mean = (t_min + t_max) / 2
    df['t_day'] = ((t_max - t_mean) * cnst.TDAY_COEF) + t_mean


def calc_prec(df: pd.DataFrame, params: dict):
    """Adjust precitation according to isoh"""
    df['prec'] *= (df.get('site_isoh', 1) / df.get('base_isoh', 1))


def calc_snowpack(df: pd.DataFrame, params: dict, snowpack=0.0):
    """Calculate snowpack as swe."""
    swe = pd.Series(snowpack, index=df.index)
    accum = (df['t_min'] <= cnst.SNOW_TCRIT)
    melt = (df['t_min'] > cnst.SNOW_TCRIT)
    swe[accum] += df['prec'][accum]
    swe[melt] -= cnst.SNOW_TRATE * (df['t_min'][melt] - cnst.SNOW_TCRIT)
    df['swe'] = np.maximum(np.cumsum(swe), 0.0)


def calc_srad_hum(df: pd.DataFrame, sg: dict, elev: float,
                  params: dict, win_type='boxcar'):
    """Calculate shortwave, humidity"""
    def _calc_tfmax(prec, dtr, sm_dtr):
        b = cnst.B0 + cnst.B1 * np.exp(-cnst.B2 * sm_dtr)
        t_fmax = 1.0 - 0.9 * np.exp(-b * np.power(dtr, cnst.C))
        inds = np.array(prec > params['sw_prec_thresh'])
        t_fmax[inds] *= cnst.RAIN_SCALAR
        return t_fmax

    # Calculate the diurnal temperature range
    df['t_max'] = np.maximum(df['t_max'], df['t_min'])
    dtr = df['t_max'] - df['t_min']
    sm_dtr = pd.Series(dtr).rolling(
        window=30, win_type=win_type, axis=0).mean().fillna(method='bfill')
    if params['n_days'] <= 30:
        warn('Timeseries is shorter than rolling mean window, filling ')
        warn('missing values with unsmoothed data')
        sm_dtr.fillna(dtr, inplace=True)

    # Effective annual prec
    if params['n_days'] <= 90:
        # Simple scaled method, minimum of 8 cm
        sum_prec = df['prec'].values.sum()
        eff_ann_prec = (sum_prec / params['n_days']) * cnst.DAYS_PER_YEAR
        eff_ann_prec = np.maximum(eff_ann_prec, 8.0)
        parray = pd.Series(eff_ann_prec, index=df.index)
    else:
        # Calculate effective annual prec using 3 month moving window
        window = pd.Series(np.zeros(params['n_days'] + 90))
        window[90:] = df['prec']

        # If yeardays at end match with those at beginning we can use
        # the end of the input to generate the beginning by looping around
        # If not, just duplicate the first 90 days
        start_day, end_day = df.index.dayofyear[0], df.index.dayofyear[-1]
        if ((start_day % 365 == (end_day % 365) + 1) or
                (start_day % 366 == (end_day % 366) + 1)):
            window[:90] = df['prec'][-90:]
        else:
            window[:90] = df['prec'][:90]

        parray = cnst.DAYS_PER_YEAR * window.rolling(
            window=90, win_type=win_type, axis=0).mean()[90:]

    # Convert to cm
    parray = parray.where(parray > 80.0, 80.0) / cnst.MM_PER_CM
    # Doing this way because parray.reindex_like(df) returns all nan
    parray.index = df.index
    df['tfmax'] = _calc_tfmax(df['prec'], dtr, sm_dtr)
    tdew = df.get('tdew', df['t_min'])
    pva = df.get('hum', svp(tdew))
    pa = atm_pres(elev)
    yday = df.index.dayofyear - 1
    df['dayl'] = sg['daylength'][yday]

    # Calculation of tdew and swrad. tdew is iterated on until
    # it converges sufficiently
    tdew_old = tdew
    tdew, pva = sw_hum_iter(df, sg, pa, pva, parray, dtr, params)
    while(np.sqrt(np.mean((tdew-tdew_old)**2)) > params['tdew_tol']):
        tdew_old = np.copy(tdew)
        tdew, pva = sw_hum_iter(df, sg, pa, pva, parray, dtr, params)
    df['vapor_pressure'] = pva


def sw_hum_iter(df, sg, pa, pva, parray, dtr, params):
    tt_max0 = sg['tt_max0']
    potrad = sg['potrad']
    daylength = sg['daylength']
    yday = df.index.dayofyear - 1

    t_tmax = np.maximum(tt_max0[yday] + (cnst.ABASE * pva), 0.0001)
    t_final = t_tmax * df['tfmax']

    # Snowpack contribution
    sc = np.zeros_like(df['swe'])
    if (cnst.MTCLIM_SWE_CORR):
        inds = np.logical_and(df['swe'] > 0.,  daylength[yday] > 0.)
        sc[inds] = ((1.32 + 0.096 * df['swe'][inds]) *
                    1.0e6 / daylength[yday][inds])
        sc = np.maximum(sc, 100.)

    # Calculation of shortwave is split into 2 components:
    # 1. Radiation from incident light
    # 2. Influence of snowpack - optionally set by MTCLIM_SWE_CORR
    df['swrad'] = potrad[yday] * t_final + sc

    # Calculate cloud effect
    if (params['lw_cloud'].upper() == 'CLOUD_DEARDORFF'):
        df['tskc'] = (1. - df['tfmax'])
    else:
        df['tskc'] = np.sqrt((1. - df['tfmax']) / 0.65)

    # Compute PET using SW radiation estimate, and update Tdew, pva **
    pet = calc_pet(df['swrad'], df['t_day'], pa, df['dayl'])
    # Calculate ratio (PET/effann_prcp) and correct the dewpoint
    ratio = pet / parray
    df['pet'] = pet * cnst.MM_PER_CM
    tmink = df['t_min'] + cnst.KELVIN
    tdew = tmink * (-0.127 + 1.121 * (1.003 - 1.444 * ratio +
                    12.312 * np.power(ratio, 2) -
                    32.766 * np.power(ratio, 3)) + 0.0006 * dtr) - cnst.KELVIN
    return tdew, svp(tdew)
