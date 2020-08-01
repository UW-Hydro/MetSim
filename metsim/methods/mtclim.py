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
from metsim.physics import atm_pres, calc_pet, svp


def run(df: pd.DataFrame, params: dict) -> pd.DataFrame:
    """
    Runs the entire mtclim set of estimation routines. This will take
    a set of daily t_min, t_max, and prec to return a set of estimated
    variables. See the documentation for the other mtclim functions for
    more information about which variables are estimated.

    Note: this function modifies the input dictionary and returns it

    Parameters
    ----------
    df:
        Dataframe containing daily inputs.
    params:
        A dictionary containing the class parameters
        of the MetSim object.

    Returns
    -------
    df:
        The same dataframe with estimated variables added
    """
    df['t_day'] = t_day(df['t_min'].values, df['t_max'].values, params)
    df['tfmax'] = tfmax(df['dtr'].values, df['smoothed_dtr'].values,
                        df['prec'].values, params)
    df['tskc'] = tskc(df['tfmax'].values, params)

    tdew_old = df['t_min'].values
    vp_temp = vapor_pressure(df['t_min'].values)
    sw_temp = shortwave(df['tfmax'].values, vp_temp,
                        df['tt_max'], df['potrad'].values)
    pet_temp = pet(sw_temp, df['t_day'].values, df['daylength'].values, params)
    tdew_temp = tdew(pet_temp, df['t_min'].values, df['seasonal_prec'].values,
                     df['dtr'].values)

    while np.sqrt(np.mean((tdew_temp - tdew_old)**2)) > params['tdew_tol']:
        tdew_old = tdew_temp.copy()
        vp_temp = vapor_pressure(tdew_temp)
        sw_temp = shortwave(df['tfmax'].values, vp_temp,
                            df['tt_max'].values, df['potrad'].values)
        pet_temp = pet(sw_temp, df['t_day'].values, df['daylength'].values,
                       params)
        tdew_temp = tdew(pet_temp, df['t_min'].values,
                         df['seasonal_prec'].values, df['dtr'].values)

    df['tdew'] = tdew_temp
    df['vapor_pressure'] = vapor_pressure(df['tdew'].values)
    df['shortwave'] = shortwave(df['tfmax'].values,
                                df['vapor_pressure'].values,
                                df['tt_max'].values, df['potrad'].values)
    df['pet'] = pet(df['shortwave'].values, df['t_day'].values,
                    df['daylength'].values, params)
    return df


def t_day(t_min: np.ndarray, t_max: np.ndarray, params: dict) -> np.ndarray:
    """
    Computes the daylight average temperature, based on a
    weighted parameterization.

    Parameters
    ----------
    t_min:
        Timeseries of daily minimum temperature
    t_max:
        Timeseries of daily maximum temperature
    params:
        Dictionary of class parameters from the MetSim object
        Note this must contain the key 'tday_coef'

    Returns
    -------
    tday:
        Daily average temperature during daylight hours
    """
    t_mean = (t_min + t_max) / 2
    return ((t_max - t_mean) * params['tday_coef']) + t_mean


def tfmax(dtr, sm_dtr, prec, params):
    """
    Computes the maximum daily transmittance of the amtosphere
    under cloudy conditions

    Parameters
    ----------
    dtr:
        Daily temperature range
    sm_dtr:
        Smoothed daily temperature range using 30 moving window
    prec:
        Daily total precipitation
    params:
        Dictionary of class parameters from the MetSim object
        Note this must contain the keys 'sw_prec_thresh' and 'rain_scalar'

    Returns
    -------
    tfmax:
        Daily maximum cloudy-sky transmittance
    """
    b = cnst.B0 + cnst.B1 * np.exp(-cnst.B2 * sm_dtr)
    tfmax = 1.0 - 0.9 * np.exp(-b * np.power(dtr, cnst.C))
    inds = np.array(prec > params['sw_prec_thresh'])
    tfmax[inds] *= params['rain_scalar']
    return tfmax


def pet(shortwave, t_day, daylength, params):
    """
    Computes potential evapotranspiration
    Note this should be jointly computed iteratively with
    ``tdew``, ``vapor_pressure``, and ``shortwave`` as used
    in the main ``run`` function.

    Parameters
    ----------
    shortwave:
        Daily estimated shortwave radiation
    t_day:
        Daylight average temperature
    daylength:
        Daily length of daylight
    params:
        Dictionary of class parameters from the MetSim object
        Note this must contain the keys 'sw_prec_thresh' and 'rain_scalar'

    Returns
    -------
    pet:
        Estimated potential evapotranspiration
    """
    pa = atm_pres(params['elev'], params['lapse_rate'])
    return calc_pet(shortwave, t_day, daylength, pa) * cnst.MM_PER_CM


def tdew(pet, t_min, seasonal_prec, dtr):
    """
    Computes dewpoint temperature
    Note this should be jointly computed iteratively with
    ``pet``, ``vapor_pressure``, and ``shortwave`` as used
    in the main ``run`` function.

    Parameters
    ----------
    pet:
        Estimated potential evapotranspiration
    t_min:
        Daily minimum temperature
    seasonal_prec:
        90 running total precipitation
    dtr:
        Daily temperature range

    Returns
    -------
    tdew:
        Estimated dewpoint temperature
    """
    parray = seasonal_prec < 80.0
    seasonal_prec[parray] = 80.0
    ratio = pet / seasonal_prec
    return np.array((t_min + cnst.KELVIN)
                    * (-0.127 + 1.121 * (1.003 - 1.444 * ratio + 12.312
                                         * np.power(ratio, 2) - 32.766
                                         * np.power(ratio, 3)) + 0.0006 * dtr)
                    - cnst.KELVIN)


def vapor_pressure(tdew):
    """
    Computes vapor pressure
    Note this should be jointly computed iteratively with
    ``pet``, ``tdew``, and ``shortwave`` as used
    in the main ``run`` function.

    Parameters
    ----------
    tdew:
        Daily dewpoint temperature

    Returns
    -------
    vapor_pressure:
        Estimated vapor pressure
    """
    return svp(tdew)


def shortwave(tfmax, vapor_pressure, tt_max, potrad):
    """
    Computes shortwave radiation
    Note this should be jointly computed iteratively with
    ``pet``, ``tdew``, and ``vapor_pressure`` as used
    in the main ``run`` function.

    Parameters
    ----------
    tfmax:
        Daily maximum cloudy-sky transmittance

    Returns
    -------
    vapor_pressure:
        Estimated vapor pressure
    """
    t_tmax = np.maximum(tt_max + (cnst.ABASE * vapor_pressure), 0.0001)
    return potrad * t_tmax * tfmax


def tskc(tfmax, params):
    """
    Computes cloud cover fraction

    Parameters
    ----------
    tfmax:
        Daily maximum cloudy-sky transmittance
    params:
        Dictionary of class parameters from the MetSim object

    Returns
    -------
    tskc:
        Daily estimated cloud cover fraction
    """
    if (params['lw_cloud'].upper() == 'CLOUD_DEARDORFF'):
        return 1. - tfmax
    return np.sqrt((1. - tfmax) / 0.65)
