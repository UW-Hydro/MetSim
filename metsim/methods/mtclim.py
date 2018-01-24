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
import xarray as xr
import xarray.ufuncs as xu

import metsim.constants as cnst
from metsim.physics import svp, calc_pet, atm_pres

# TODO:
# - Why isn't xu.power a thing
# - sg needs to be a xarray dataset that is broadcastable with ds
# - use scipy.optimize instead of homegrown while loop to iterate when solving
#   for t_dew
# - fix namespace conflict for calc_pet


def run(ds: xr.Dataset, params: dict, sg: xr.Dataset):
    '''
    Run the MTCLIM method

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with t_min, t_max, dtr, smoothed_dtr, and prec in it
    params : dict
        Parameters for MTCLIM
    sg : xr.Dataset
        Solar geometry dataset.

    Returns
    -------
    ds : xr.Dataset
        Dataset with same variables as ds, plus tfmax, tskc, dayl, tdew,
        vapor_pressure, shortwave, and pet.
    '''

    ds['t_day'] = calc_t_day(ds['t_min'], ds['t_max'], params)
    ds['tfmax'] = calc_tfmax(ds['dtr'], ds['smoothed_dtr'], ds['prec'], params)
    ds['tskc'] = calc_tskc(ds['tfmax'], params)
    yday = ds.indexes['time'].dayofyear - 1
    ds['dayl'] = calc_dayl(yday, sg)

    vp_temp = svp(ds['t_min'])
    sw_temp = calc_shortwave(yday, ds['tfmax'], vp_temp, sg)
    pet_temp = _calc_pet(sw_temp, ds['t_day'], ds['dayl'], params)
    tdew_temp = calc_tdew(pet_temp, ds['t_min'], ds['seasonal_prec'],
                          ds['dtr'])

    # Adjust seasonal precipitation
    seasonal_prec = ds['seasonal_prec'].where(ds['seasonal_prec'] >= 80., 80)

    ds['tdew'] = calc_tdew_iter(tdew_temp, yday, ds['tfmax'],
                                ds['t_day'], ds['dayl'], ds['t_min'],
                                seasonal_prec, ds['dtr'], params, sg)
    ds['vapor_pressure'] = svp(ds['tdew'])
    ds['shortwave'] = calc_shortwave(
        yday, ds['tfmax'], ds['vapor_pressure'], sg)
    ds['pet'] = _calc_pet(ds['shortwave'], ds['t_day'], ds['dayl'], params)
    return ds


def _calc_tdew_iter(tdew_temp, yday, tfmax, t_day, dayl, t_min,
                    seasonal_prec, dtr, params, sg):
    '''helper function to facilitate calculating tdew

    this function must be vectorized and applied to each timestep
    '''

    tdew_old = t_min

    while(np.sqrt(np.mean((tdew_temp-tdew_old)**2)) > params['tdew_tol']):
        tdew_old = tdew_temp.copy()
        vp = svp(tdew_temp)
        sw = calc_shortwave(yday, tfmax, vp, sg)
        pet = _calc_pet(sw, t_day, dayl, params)
        tdew_temp = calc_tdew(pet, t_min, seasonal_prec, dtr)

    return tdew_temp


def calc_tdew_iter(*args, **kwargs):
    return xr.apply_ufunc(_calc_tdew_iter, *args, **kwargs,
                          dask='parallized', vectorized=True)


def calc_t_day(t_min, t_max, params):
    t_mean = (t_min + t_max) / 2
    return ((t_max - t_mean) * params['tday_coef']) + t_mean


def calc_tfmax(dtr, sm_dtr, prec, params):
    b = cnst.B0 + cnst.B1 * xu.exp(-cnst.B2 * sm_dtr)
    tfmax = 1.0 - 0.9 * xu.exp(-b * np.power(dtr, cnst.C))
    tfmax = tfmax.where(prec <= params['sw_prec_thresh'],
                        tfmax * params['rain_scalar'])
    return tfmax


def calc_dayl(yday, sg):
    return sg['daylength'][yday]


def _calc_pet(shortwave, t_day, dayl, params):
    pa = atm_pres(params['elev'], params['lapse_rate'])
    return calc_pet(shortwave, t_day, dayl, pa) * cnst.MM_PER_CM


def calc_tdew(pet, t_min, seasonal_prec, dtr):
    ratio = pet / seasonal_prec
    return ((t_min + cnst.KELVIN) * (-0.127 + 1.121 * (1.003 - 1.444 * ratio
            + 12.312 * np.power(ratio, 2) - 32.766 * np.power(ratio, 3))
            + 0.0006 * dtr) - cnst.KELVIN)


def calc_shortwave(yday, tfmax, vapor_pressure, sg):
    t_tmax = xu.maximum(sg['tt_max0'][yday]
                        + (cnst.ABASE * vapor_pressure), 0.0001)
    return sg['potrad'][yday] * t_tmax * tfmax


def calc_tskc(tfmax, params):
    if (params['lw_cloud'].upper() == 'CLOUD_DEARDORFF'):
        return 1. - tfmax
    return xu.sqrt((1. - tfmax) / 0.65)
