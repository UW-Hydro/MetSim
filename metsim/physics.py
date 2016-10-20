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
from metsim.defaults import CONSTS as constants


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
    gamma = constants['CP'] * pa / (lhvap * constants['EPS'])

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
    t1 = 1.0 - (constants['LR_STD'] * elev) / constants['T_STD']
    t2 = constants['G_STD'] / (constants['LR_STD'] * (constants['R'] /
                                                      constants['MA']))

    return constants['P_STD'] * np.power(t1, t2)


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
