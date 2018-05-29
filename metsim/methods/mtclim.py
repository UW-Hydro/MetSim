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

import metsim.constants as cnst
from metsim.physics import atm_pres, calc_pet, svp


def run(df, params):
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

    while(np.sqrt(np.mean((tdew_temp - tdew_old)**2)) > params['tdew_tol']):
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


def t_day(t_min, t_max, params):
    t_mean = (t_min + t_max) / 2
    return ((t_max - t_mean) * params['tday_coef']) + t_mean


def tfmax(dtr, sm_dtr, prec, params):
    b = cnst.B0 + cnst.B1 * np.exp(-cnst.B2 * sm_dtr)
    tfmax = 1.0 - 0.9 * np.exp(-b * np.power(dtr, cnst.C))
    inds = np.array(prec > params['sw_prec_thresh'])
    tfmax[inds] *= params['rain_scalar']
    return tfmax


def pet(shortwave, t_day, daylength, params):
    pa = atm_pres(params['elev'], params['lapse_rate'])
    return calc_pet(shortwave, t_day, daylength, pa) * cnst.MM_PER_CM


def tdew(pet, t_min, seasonal_prec, dtr):
    parray = seasonal_prec < 80.0
    seasonal_prec[parray] = 80.0
    ratio = pet / seasonal_prec
    return ((t_min + cnst.KELVIN) * (-0.127 + 1.121 * (1.003 - 1.444 * ratio + 12.312 * \
            np.power(ratio, 2) - 32.766 * np.power(ratio, 3)) + 0.0006 * dtr) - cnst.KELVIN)


def vapor_pressure(tdew):
    return svp(tdew)


def shortwave(tfmax, vapor_pressure, tt_max, potrad):
    t_tmax = np.maximum(tt_max + (cnst.ABASE * vapor_pressure), 0.0001)
    return potrad * t_tmax * tfmax


def tskc(tfmax, params):
    if (params['lw_cloud'].upper() == 'CLOUD_DEARDORFF'):
        return 1. - tfmax
    return np.sqrt((1. - tfmax) / 0.65)
