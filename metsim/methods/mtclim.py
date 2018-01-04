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


def run(df: pd.DataFrame, params: dict, sg: dict):
    """
    Run all of the mtclim df generation

    Parameters
    ----------
    df: pd.DataFrame
        The daily dfs given from input
    params: dict
        Dictionary of parameters from a MetSim object
    solar_geom: dict
        Solar geometry of the site

    Returns
    -------
    df:
        Dataframe of daily or subdaily dfs
    """
    t_mean = (df['t_min'] + df['t_max']) / 2
    df['t_day'] = (
            (df['t_max'] - t_mean) * params['tday_coef']) + t_mean
    b = cnst.B0 + cnst.B1 * np.exp(-cnst.B2 * df['smoothed_dtr'])
    df['tfmax'] = (
            1.0 - 0.9 * np.exp(-b * np.power(df['dtr'], cnst.C)))
    inds = np.array(df['prec'] > params['sw_prec_thresh'])
    df['tfmax'][inds] *= params['rain_scalar']

    tdew = df.get('tdew', df['t_min'])
    df['vapor_pressure'] = df.get('vapor_pressure', svp(tdew.values))
    pa = atm_pres(params['elev'], params['lapse_rate'])
    yday = df.index.dayofyear - 1
    df['dayl'] = sg['daylength'][yday]

    # Calculation of tdew and shortwave. tdew is iterated on until
    # it converges sufficiently
    tdew_old = tdew
    err = 500
    while(err > params['tdew_tol']):
        tdew_old = np.copy(tdew)
        yday = df.index.dayofyear - 1
        t_tmax = np.maximum(sg['tt_max0'][yday]
                            + (cnst.ABASE * df['vapor_pressure']), 0.0001)
        df['shortwave'] = sg['potrad'][yday] * t_tmax * df['tfmax']
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
        err = np.sqrt(np.mean((tdew-tdew_old)**2))
    df['vapor_pressure'] = svp(tdew.values)

    # Calculate cloud effect
    if (params['lw_cloud'].upper() == 'CLOUD_DEARDORFF'):
        df['tskc'] = (1. - df['tfmax'])
    else:
        df['tskc'] = np.sqrt((1. - df['tfmax']) / 0.65)
    return df
