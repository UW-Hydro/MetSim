"""
TODO
"""

import numpy as np
import pandas as pd 

from metsim.configuration import PARAMS as params
from metsim.configuration import CONSTS as consts
from metsim.configuration import OPTIONS as options

class Forcing(object):
    """
    TODO
    """
    
    def __init__(self, met_data: pd.DataFrame, params: dict):
        """
        Initialize the forcing object.

        Args:
            met_data: A dataframe containing given forcings.
                      For example, these could be read in 
                      from a file or generated by some other
                      code.
            params: A dictionary of parameters about the site
                    being analyzed.
        """
        # Update the configuration
        import metsim.configuration
        metsim.configuration.update(params)
        self._validate_inputs(met_data, params)
        self.met_data = met_data
        self.solar_geom = self._calc_solar_geom()
        

    def generate_met_forcings(self, forcing_method):
        """
        Generate the required forcings.  In the future this 
        should have some added functionality so that it can
        pick out the required forcings if there are more than
        the minimum requirement given.  
        """
        # TODO: This is a shortcut, but works for the default
        #       situation
        forcing_method.run(self.met_data, self.solar_geom)


    def disaggregate(self):
        """
        Converts all of the daily variables to hourly
        """
        import metsim.disaggregate as disagg
        self.met_data = disagg.disaggregate(self.met_data, self.solar_geom) 


    def set_dates(self, dates: pd.DatetimeIndex):
        """ Set the index of the met_data to certain dates """
        self.met_data.set_index(dates, inplace=True)
        self.add_data("day_of_year", dates.dayofyear)


    def add_data(self, key, data):
        """ Add some data to the met_data with a given key """
        self.met_data[key] = data


    def _validate_inputs(self, met_data: pd.DataFrame, params: dict):
        """ Make sure all the required information has been input """
        met_data['wind']
        met_data['precip']
        met_data['t_min']
        met_data['t_max']
        params['lat']
        params['lon']


    def _calc_solar_geom(self):
        """
        Flat earth assumption
        """
        dayperyear = np.ceil(consts['DAYS_PER_YEAR'])
        tt_max0   = np.zeros(dayperyear)
        daylength = np.zeros(dayperyear)
        potrad    = np.zeros(dayperyear) 
        t1 = 1.0 - (consts['LR_STD'] * params['site_elev'])/consts['T_STD']
        t1 = 1.0 - (consts['LR_STD'] * 1668.15)/consts['T_STD']
        t2 = consts['G_STD'] / (consts['LR_STD'] * (consts['R'] / consts['MA']))
        trans = np.power(params['TBASE'], np.power(t1, t2))
        print(params['site_elev'])
        
        lat    = np.clip(params['lat']*consts['RADPERDEG'], -np.pi/2., np.pi/2.0)
        coslat = np.cos(lat)
        sinlat = np.sin(lat)
        dt     = consts['SRADDT']  
        dh     = dt / consts['SECPERRAD'] 
    
        tiny_step_per_day = consts['SEC_PER_DAY'] / consts['SRADDT']
        tiny_rad_fract    = np.zeros(shape=(dayperyear, tiny_step_per_day), dtype=np.float64)
        for i in range(int(dayperyear-1)):
            decl = consts['MINDECL'] * np.cos((i + consts['DAYSOFF']) * consts['RADPERDAY'])
            cosdecl = np.cos(decl)
            sindecl = np.sin(decl)
            cosegeom = coslat * cosdecl

            sinegeom = sinlat * sindecl
            coshss = np.clip(-sinegeom / cosegeom, -1, 1)
            hss = np.arccos(coshss)  
            daylength[i] = np.minimum(2.0 * hss * consts['SECPERRAD'], consts['SEC_PER_DAY'])
            dir_beam_topa = (1368.0 + 45.5 * np.sin((2.0 * np.pi * i / 365.25) + 1.7)) * dt
            h = np.arange(-hss, hss, dh)
            cosh = np.cos(h)
            cza  = cosegeom * cosh + sinegeom
            cza_inds = np.array(cza > 0.)
           
            dir_flat_topa = -1 * np.ones(len(h))
            dir_flat_topa[cza_inds] = dir_beam_topa * cza[cza_inds]
            
            am = np.zeros(len(h))
            am[cza_inds] = 1.0 / (cza[cza_inds] + 0.0000001)
            am_inds = np.array(am > 2.9)

            ami = np.zeros(len(am_inds))
            ami = (np.arccos(cza[am_inds])/consts['RADPERDEG'] - 69).astype(int)
            if len(ami) != 0:
                ami = np.clip(ami, 0, 20)
                am[am_inds] = consts['OPTAM'][ami]
    
            trans2 = np.power(trans, am)
            
            sum_trans = sum(trans2 * dir_flat_topa)
            sum_flat_potrad = sum(dir_flat_topa)
            tinystep = np.clip((12*consts['SEC_PER_HOUR']+h*consts['SECPERRAD'])/dt,
                            0, tiny_step_per_day-1).astype(int)
            tiny_rad_fract[i][tinystep] = dir_flat_topa
            
            if daylength[i] and sum_flat_potrad > 0:
                tiny_rad_fract[i] /= sum_flat_potrad

            if daylength[i] > 0:
                tt_max0[i] = sum_trans / sum_flat_potrad
                potrad[i] = sum_flat_potrad / daylength[i]
            else:
                tt_max0[i] = 0.
                potrad[i] = 0.
        tt_max0[-1] = tt_max0[-2]
        potrad[-1] = potrad[-2]
        daylength[-1] = daylength[-2]
        tiny_rad_fract[-1] = tiny_rad_fract[-2]
        solar_geom = {"tt_max0" : tt_max0,
                      "potrad" : potrad,
                      "daylength" : daylength,
                      "tiny_rad_fract" : tiny_rad_fract}
        return solar_geom 
    
