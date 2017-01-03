'''
physics
'''

# Mountain Climate Simulator, meteorological forcing disaggregator
# Copyright (C) 2015  Joe Hamman

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
from metsim.configuration import CONSTS as consts
from metsim.configuration import PARAMS as params

def calc_pet(rad, ta, pa, dayl, dt=0.2):
    '''
    calculates the potential evapotranspiration for aridity corrections in
    `calc_vpd()`, according to Kimball et al., 1997

    Parameters
    ----------
    rad : scalar or numpy.ndarray
        daylight average incident shortwave radiation (W/m2)
    ta : scalar or numpy.ndarray
        daylight average air temperature (deg C)
    pa : scalar or numpy.ndarray
        air pressure (Pa)
    dayl : scalar or numpy.ndarray
        daylength (s)
    dt : scalar, optional
        offset for saturation vapor pressure calculation, default = 0.2

    Returns
    ----------
    pet : scalar or numpy.ndarray
        Potential evapotranspiration (cm/day)
    '''
    # rnet       # (W m-2) absorbed shortwave radiation avail. for ET
    # lhvap      # (J kg-1) latent heat of vaporization of water
    # gamma      # (Pa K-1) psychrometer parameter
    # dt = 0.2   # offset for saturation vapor pressure calculation
    # t1, t2     # (deg C) air temperatures
    # pvs1, pvs2 # (Pa)   saturated vapor pressures
    # pet        # (kg m-2 day-1) potential evapotranspiration
    # s          # (Pa K-1) slope of saturated vapor pressure curve

    # calculate absorbed radiation, assuming albedo = 0.2  and ground
    # heat flux = 10% of absorbed radiation during daylight
    rnet = rad * 0.72

    # calculate latent heat of vaporization as a function of ta
    lhvap = 2.5023e6 - 2430.54 * ta

    # calculate the psychrometer parameter: gamma = (cp pa)/(lhvap epsilon)
    # where:
    # cp       (J/kg K)   specific heat of air
    # epsilon  (unitless) ratio of molecular weights of water and air
    gamma = consts['CP'] * pa / (lhvap * consts['EPS'])

    # estimate the slope of the saturation vapor pressure curve at ta
    # temperature offsets for slope estimate
    t1 = ta + dt
    t2 = ta - dt

    # calculate saturation vapor pressures at t1 and t2, using formula from
    # Abbott, P.F., and R.C. Tabony, 1985. The estimation of humidity
    # parameters. Meteorol. Mag., 114:49-56.
    pvs1 = svp(t1)
    pvs2 = svp(t2)

    # calculate slope of pvs vs. T curve near ta
    s = (pvs1 - pvs2) / (t1 - t2)
    # can this be s = svp_slope(ta)? JJH

    # calculate PET using Priestly-Taylor approximation, with coefficient
    # set at 1.26. Units of result are kg/m^2/day, equivalent to mm water/day
    pet = (1.26 * (s / (s + gamma)) * rnet * dayl) / lhvap

    # return a value in centimeters/day, because this value is used in a ratio
    # to annual total precip, and precip units are centimeters
    return (pet / 10.)


def atm_pres(elev):
    '''atmospheric pressure (Pa) as a function of elevation (m)

    Parameters
    ----------
    elev : scalar or numpy.ndarray
        Elevation (meters)

    Returns
    -------
    pressure : scalar or numpy.ndarray
        Atmospheric pressure at elevation `elev` (Pa)

    References
    ----------
    * Iribane, J.V., and W.L. Godson, 1981. Atmospheric Thermodynamics, 2nd
      Edition. D. Reidel Publishing Company, Dordrecht, The Netherlands.
      (p. 168)
    '''
    t1 = 1.0 - (consts['LR_STD'] * elev) / consts['T_STD']
    t2 = consts['G_STD'] / (consts['LR_STD'] * (consts['R']/consts['MA']))
    return consts['P_STD'] * np.power(t1, t2)


def svp(temp, a=0.61078, b=17.269, c=237.3):
    '''Compute the saturated vapor pressure.

    Parameters
    ----------
    temp : numpy.ndarray
        Temperature (degrees Celsius)

    Returns
    ----------
    pressure : numpy.ndarray
        Saturated vapor pressure at temperature `temp` (Pa)

    References
    ----------
    * Maidment, David R. Handbook of hydrology. McGraw-Hill Inc., 1992.
      Equation 4.2.2.
    '''
    svp = a * np.exp((b * temp) / (c + temp))
    inds = np.nonzero(temp < 0.)[0]
    svp[inds] *= 1.0 + .00972 * temp[inds] + .000042 * np.power(temp[inds], 2)
    return svp * 1000.


def svp_slope(temp, a=0.61078, b=17.269, c=237.3):
    '''Compute the gradient of the saturated vapor pressure as a function of
    temperature.

    Parameters
    ----------
    temp : numpy.ndarray
        Temperature (degrees Celsius)

    Returns
    -------
    gradient : numpy.ndarray
        Gradient of d(svp)/dT.

    References
    ----------
    * Maidment, David R. Handbook of hydrology. McGraw-Hill Inc., 1992.
      Equation 4.2.3.
    '''
    return (b * c) / ((c + temp) * (c + temp)) * svp(temp, a=a, b=b, c=c)


def solar_geom(elev: float, lat: float) -> dict:
    """
    Flat earth assumption
    """
    dayperyear = int(np.ceil(consts['DAYS_PER_YEAR']))
    tt_max0   = np.zeros(dayperyear)
    daylength = np.zeros(dayperyear)
    potrad    = np.zeros(dayperyear) 
    t1 = 1.0 - (consts['LR_STD'] * elev)/consts['T_STD']
    t2 = consts['G_STD'] / (consts['LR_STD'] * (consts['R'] / consts['MA']))
    trans = np.power(params['TBASE'], np.power(t1, t2))
   
    # Translate lat to rad
    lat    = np.clip(lat * consts['RADPERDEG'], -np.pi/2., np.pi/2.0)
    coslat = np.cos(lat)
    sinlat = np.sin(lat)

    # Sub-daily time step and angular step
    dt     = consts['SRADDT']  
    dh     = dt / consts['SECPERRAD'] 

    tiny_step_per_day = int(consts['SEC_PER_DAY'] / consts['SRADDT'])
    tiny_rad_fract    = np.zeros(shape=(dayperyear, tiny_step_per_day), dtype=np.float64)
    for i in range(int(dayperyear-1)):
        # Declination and quantities of interest
        decl = consts['MINDECL'] * np.cos((i + consts['DAYSOFF']) * consts['RADPERDAY'])
        cosdecl = np.cos(decl)
        sindecl = np.sin(decl)
        cosegeom = coslat * cosdecl
    
        sinegeom = sinlat * sindecl
        coshss = np.clip(-sinegeom / cosegeom, -1, 1)
        hss = np.arccos(coshss)  
        daylength[i] = np.minimum(2.0 * hss * consts['SECPERRAD'], consts['SEC_PER_DAY'])
        dir_beam_topa = (1368.0 + 45.5 * np.sin((2.0 * np.pi * i / 365.25) + 1.7)) * dt

        # Set up angular calculations in vectorized array
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
 
