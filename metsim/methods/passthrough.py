"""
Passthrough
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
from metsim.methods.mtclim import t_day, tfmax, tskc, pet, tdew, vapor_pressure

def run(df, params):
    assert 'shortwave' in df

    if 't_day' not in df:
        df['t_day'] = t_day(df['t_min'].values, df['t_max'].values, params)
    if 'tfmax' not in df:
        df['tfmax'] = tfmax(df['dtr'].values, df['smoothed_dtr'].values,
                            df['prec'].values, params)
    if 'tskc' not in df:
        df['tskc'] = tskc(df['tfmax'].values, params)
    if 'pet' not in df:
        df['pet'] = pet(df['shortwave'].values, df['t_day'].values,
                        df['daylength'].values, params)
    if 'tdew' not in df:
        df['tdew'] = tdew(df['pet'].values, df['t_min'].values,
                          df['seasonal_prec'].values, df['dtr'].values)
    if 'vapor_pressure' not in df:
        df['vapor_pressure'] = vapor_pressure(df['tdew'].values)
    return df
