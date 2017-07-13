"""
MTCLIM
"""
# Meteorology Simulator
# Copyright (C) 2017  The Computational Hydrology Group, Department of Civil
# and Environmental Engineering, University of Washington.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import pandas as pd

import metsim.constants as cnst
from metsim.physics import svp, calc_pet, atm_pres


def run(forcing: pd.DataFrame, params: dict, sg: dict,
        elev: float, swe: float):
    """
    Run all of the mtclim forcing generation

    Parameters
    ----------
    forcing: pd.DataFrame
        The daily forcings given from input
    params: dict
        Dictionary of parameters from a MetSim object
    solar_geom: dict
        Solar geometry of the site
    elev: float
        Elevation of the site being simulated
    swe: float
        Initial snow water equivalent (SWE) for
        the site being simulated

    Returns
    -------
    forcing:
        Dataframe of daily or subdaily forcings
    """
    params['n_days'] = len(forcing)
    calc_t_air(forcing, elev, params)
    calc_snowpack(forcing, params, swe)
    calc_srad_hum(forcing, sg, elev, params)

    return forcing


def calc_t_air(df: pd.DataFrame, elev: float, params: dict):
    """
    Calculate mean daily temperature

    Parameters
    ----------
    df:
        Dataframe with daily max and min temperatures
    elev:
        Elevation in meters
    params:
        Dictionary containing parameters from a
        MetSim object.
    """
    t_mean = (df['t_min'] + df['t_max']) / 2
    df['t_day'] = ((df['t_max'] - t_mean) * params['tday_coef']) + t_mean


def calc_snowpack(df: pd.DataFrame, params: dict, snowpack: float=0.0):
    """
    Estimate snowpack as swe.

    Parameters
    ----------
    df:
        Dataframe with daily timeseries of precipitation
        and minimum temperature.
    snowpack:
        (Optional - defaults to 0) Initial snowpack
    """
    swe = pd.Series(snowpack, index=df.index)
    accum = (df['t_min'] <= params['snow_crit_temp'])
    melt = (df['t_min'] > params['snow_crit_temp'])
    swe[accum] += df['prec'][accum]
    swe[melt] -= (params['snow_melt_rate']
                  * (df['t_min'][melt] - params['snow_crit_temp']))
    df['swe'] = np.maximum(np.cumsum(swe), 0.0)


def calc_srad_hum(df: pd.DataFrame, sg: dict, elev: float,
                  params: dict):
    """
    Calculate shortwave, humidity

    Parameters
    ----------
    df:
        Dataframe containing daily timeseries
    elev:
        Elevation in meters
    params:
        A dictionary of parameters from the
        MetSim object
    """
    def _calc_tfmax(prec, dtr, sm_dtr):
        """Estimate cloudy day transmittance"""
        b = cnst.B0 + cnst.B1 * np.exp(-cnst.B2 * sm_dtr)
        t_fmax = 1.0 - 0.9 * np.exp(-b * np.power(dtr, cnst.C))
        inds = np.array(prec > params['sw_prec_thresh'])
        t_fmax[inds] *= params['rain_scalar']
        return t_fmax

    # Calculate the diurnal temperature range
    df['t_max'] = np.maximum(df['t_max'], df['t_min'])
    dtr = df['t_max'] - df['t_min']
    df['dtr'] = dtr
    sm_dtr = df['smoothed_dtr']
    df['tfmax'] = _calc_tfmax(df['prec'], dtr, sm_dtr)
    tdew = df.get('tdew', df['t_min'])
    pva = df.get('hum', svp(tdew.values))
    pa = atm_pres(elev, params['lapse_rate'])
    yday = df.index.dayofyear - 1
    df['dayl'] = sg['daylength'][yday]

    # Calculation of tdew and shortwave. tdew is iterated on until
    # it converges sufficiently
    tdew_old = tdew
    tdew, pva = sw_hum_iter(df, sg, pa, pva, dtr, params)
    while(np.sqrt(np.mean((tdew-tdew_old)**2)) > params['tdew_tol']):
        tdew_old = np.copy(tdew)
        tdew, pva = sw_hum_iter(df, sg, pa, pva, dtr, params)
    df['vapor_pressure'] = pva


def sw_hum_iter(df: pd.DataFrame, sg: dict, pa: float, pva: pd.Series,
                dtr: pd.Series, params: dict):
    """
    Calculated updated values for dewpoint temperature
    and saturation vapor pressure.

    Parameters
    ----------
    df:
        Dataframe containing daily timeseries of
        cloud cover fraction, tfmax, swe, and
        shortwave radiation
    sg:
        Solar geometry dictionary, calculated with
        `metsim.physics.solar_geom`.
    pa:
        Air pressure in Pascals
    pva:
        Vapor presure in Pascals
    dtr:
        Daily temperature range
    params:
        A dictionary of parameters from a MetSim object

    Returns
    -------
    (tdew, svp):
        A tuple of dewpoint temperature and saturation
        vapor pressure
    """
    tt_max0 = sg['tt_max0']
    potrad = sg['potrad']
    daylength = sg['daylength']
    yday = df.index.dayofyear - 1

    t_tmax = np.maximum(tt_max0[yday] + (cnst.ABASE * pva), 0.0001)
    t_final = t_tmax * df['tfmax']
    df['pva'] = pva
    df['tfinal'] = t_final

    # Snowpack contribution
    sc = np.zeros_like(df['swe'])
    if (params['mtclim_swe_corr']):
        inds = np.logical_and(df['swe'] > 0.,  daylength[yday] > 0.)
        sc[inds] = ((1.32 + 0.096 * df['swe'][inds]) *
                    1.0e6 / daylength[yday][inds])
        sc = np.maximum(sc, 100.)

    # Calculation of shortwave is split into 2 components:
    # 1. Radiation from incident light
    # 2. Influence of snowpack - optionally set by MTCLIM_SWE_CORR
    df['shortwave'] = potrad[yday] * t_final + sc

    # Calculate cloud effect
    if (params['lw_cloud'].upper() == 'CLOUD_DEARDORFF'):
        df['tskc'] = (1. - df['tfmax'])
    else:
        df['tskc'] = np.sqrt((1. - df['tfmax']) / 0.65)

    # Compute PET using SW radiation estimate, and update Tdew, pva **
    pet = calc_pet(df['shortwave'].values, df['t_day'].values,
                   df['dayl'].values, pa)
    # Calculate ratio (PET/effann_prcp) and correct the dewpoint
    parray = df['seasonal_prec'] / cnst.MM_PER_CM
    ratio = pet / parray.where(parray > 8.0, 8.0)
    df['pet'] = pet * cnst.MM_PER_CM
    tmink = df['t_min'] + cnst.KELVIN
    tdew = tmink * (-0.127 + 1.121 * (1.003 - 1.444 * ratio +
                    12.312 * np.power(ratio, 2) -
                    32.766 * np.power(ratio, 3)) + 0.0006 * dtr) - cnst.KELVIN
    return tdew, svp(tdew.values)
