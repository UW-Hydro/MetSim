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


def run(forcing: pd.DataFrame, params: dict, sg: dict):
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

    Returns
    -------
    forcing:
        Dataframe of daily or subdaily forcings
    """
    params['n_days'] = len(forcing)
    t_mean = (forcing['t_min'] + forcing['t_max']) / 2
    forcing['t_day'] = (
            (forcing['t_max'] - t_mean) * params['tday_coef']) + t_mean
    b = cnst.B0 + cnst.B1 * np.exp(-cnst.B2 * forcing['smoothed_dtr'])
    forcing['tfmax'] = (
            1.0 - 0.9 * np.exp(-b * np.power(forcing['dtr'], cnst.C)))
    inds = np.array(forcing['prec'] > params['sw_prec_thresh'])
    forcing['tfmax'][inds] *= params['rain_scalar']

    tdew = forcing.get('tdew', forcing['t_min'])
    forcing['vapor_pressure'] = forcing.get('vapor_pressure', svp(tdew.values))
    pa = atm_pres(params['elev'], params['lapse_rate'])
    yday = forcing.index.dayofyear - 1
    forcing['dayl'] = sg['daylength'][yday]

    # Calculation of tdew and shortwave. tdew is iterated on until
    # it converges sufficiently
    tdew_old = tdew
    tdew = sw_hum_iter(forcing, sg, pa, params)
    while(np.sqrt(np.mean((tdew-tdew_old)**2)) > params['tdew_tol']):
        tdew_old = np.copy(tdew)
        tdew = sw_hum_iter(forcing, sg, pa, params)
    forcing['vapor_pressure'] = svp(tdew.values)
    return forcing


def sw_hum_iter(df: pd.DataFrame, sg: dict, pa: float, params: dict):
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
    params:
        A dictionary of parameters from a MetSim object

    Returns
    -------
    Dewpoint temperature
    """
    yday = df.index.dayofyear - 1
    t_tmax = np.maximum(sg['tt_max0'][yday]
                        + (cnst.ABASE * df['vapor_pressure']), 0.0001)
    df['shortwave'] = sg['potrad'][yday] * t_tmax * df['tfmax']

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
    tdew = (tmink * (-0.127 + 1.121 * (1.003 - 1.444 * ratio
            + 12.312 * np.power(ratio, 2) - 32.766 * np.power(ratio, 3))
            + 0.0006 * df['dtr']) - cnst.KELVIN)
    return tdew
