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
from numba import jit
import metsim.configuration as conf

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
    gamma = conf.CP * pa / (lhvap * conf.EPS)

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
    t1 = 1.0 - (conf.LR_STD * elev) / conf.T_STD
    t2 = conf.G_STD / (conf.LR_STD * (conf.R/conf.MA))
    return conf.P_STD * np.power(t1, t2)


@jit
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




@jit(nopython=True, cache=True)
def solar_geom(elev, lat):
    """
    Flat earth assumption
    """
    # optical airmass by degrees
    OPTAM = [2.90,  3.05,  3.21,  3.39,  3.69,  3.82,  4.07, 
             4.37,  4.72,  5.12,  5.60,  6.18,  6.88,  7.77, 
             8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00]
    dayperyear = int(np.ceil(conf.DAYS_PER_YEAR))
    tt_max0   = np.zeros(dayperyear)
    daylength = np.zeros(dayperyear)
    flat_potrad    = np.zeros(dayperyear) 
    t1 = 1.0 - (conf.LR_STD * elev)/conf.T_STD
    t2 = conf.G_STD / (conf.LR_STD * (conf.R /conf.MA))
    trans = np.power(conf.TBASE, np.power(t1, t2))
   
    # Translate lat to rad
    lat    = np.minimum(np.maximum(lat * conf.RADPERDEG, -np.pi/2.), np.pi/2.0)
    coslat = np.cos(lat)
    sinlat = np.sin(lat)

    # Sub-daily time step and angular step
    dt = conf.SRADDT  
    dh = dt / conf.SECPERRAD 

    tiny_step_per_day = int(conf.SEC_PER_DAY / conf.SRADDT)
    tiny_rad_fract = np.zeros((dayperyear, tiny_step_per_day))
    for i in range(dayperyear-1):
        # Declination and quantities of interest
        decl = conf.MINDECL * np.cos((i + conf.DAYSOFF) * conf.RADPERDAY)
        cosdecl = np.cos(decl)
        sindecl = np.sin(decl)
       
        # calculate daylength as a function of lat and decl
        cosegeom = coslat * cosdecl
        sinegeom = sinlat * sindecl
        coshss = min(max(-sinegeom / cosegeom, -1), 1)
        hss = np.arccos(coshss)  
        daylength[i] = min(2.0 * hss * conf.SECPERRAD, conf.SEC_PER_DAY)
        dir_beam_topa = (1368.0+45.5*np.sin((2.0*np.pi*i/conf.DAYS_PER_YEAR)+1.7))*dt
        sum_trans = 0
        sum_flat_potrad= 0
        # Set up angular calculations 
        for h in np.arange(-hss, hss, dh):
            cosh = np.cos(h)
            cza  = cosegeom * cosh + sinegeom
            if (cza > 0):
                dir_flat_topa = dir_beam_topa * cza
                am = 1.0 / (cza + 0.0000001)
                if (am > 2.9):
                    ami = min(max(int(np.arccos(cza)/conf.RADPERDEG) - 69, 0), 20)
                    am = OPTAM[ami]
                sum_trans += (np.power(trans, am) * dir_flat_topa)
                sum_flat_potrad += dir_flat_topa
            else:
                dir_flat_topa = 0 
            
            tinystep = int(min(
                               max((12*conf.SEC_PER_HOUR+h*conf.SECPERRAD)/dt,0), 
                               tiny_step_per_day-1))
            tiny_rad_fract[i][tinystep] = dir_flat_topa
            
        if daylength[i] and sum_flat_potrad > 0:
           tiny_rad_fract[i] /= sum_flat_potrad

        if daylength[i]:
            tt_max0[i] = sum_trans / sum_flat_potrad
            flat_potrad[i] = sum_flat_potrad / daylength[i]
        else:
            tt_max0[i] = 0.
            flat_potrad[i] = 0.
    tt_max0[dayperyear-1] = tt_max0[dayperyear-2]
    flat_potrad[dayperyear-1] = flat_potrad[dayperyear-2]
    daylength[dayperyear-1] = daylength[dayperyear-2]
    tiny_rad_fract[dayperyear-1] = tiny_rad_fract[dayperyear-2]

    return tiny_rad_fract, daylength, flat_potrad, tt_max0


