"""
Stores default constants, parameters, and options.  
May be overwritten from config files.
"""
import numpy as np

CONSTS = {
          'SECPERRAD': 13750.9871,  # seconds per radian of hour angle
          'RADPERDAY': 0.017214,    # radians of Earth orbit per julian day
          'RADPERDEG': 0.01745329,  # radians per degree
          'MINDECL': -0.4092797,    # minimum declination (radians)
          'DAYSOFF': 11.25,         # julian day offset of winter solstice
          'SRADDT': 30.0,           # timestep for radiation routine (seconds)
                                    # Note:  Make sure that 3600 % SRADDT == 0
          'MA': 28.9644e-3,         # (kg mol-1) molecular weight of air
          'R': 8.3143,              # (m3 Pa mol-1 K-1) gas law constant
          'G_STD': 9.80665,         # (m s-2) standard gravitational accel.
          'P_STD': 101325.0,        # (Pa) standard pressure at 0. m elevation
          'T_STD': 288.15,          # (K) standard temp at 0. m elevation
          'CP': 1010.0,             # (J kg-1 K-1) specific heat of air
          'LR_STD': 0.0065,         # (-K m-1) standard temperature lapse rate
           
          # optical airmass by degrees
          'OPTAM': np.array([2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07, 4.37,
                            4.72, 5.12, 5.60, 6.18, 6.88, 7.77, 8.90, 10.39,
                            12.44, 15.36, 19.79, 26.96, 30.00]),
          'KELVIN': 273.15,
          'EPS': 0.62196351,
          'DAYS_PER_YEAR' : 365.25,
          'SEC_PER_DAY'   : 86400,
          'HOURS_PER_DAY'  : 24,
          'SEC_PER_HOUR'  : 3600,
          'STEFAN_B' : 5.669e-8  # (W m^-22 K^-4) Stefan Boltzmann constant 
        }

PARAMS = {
          'TDAYCOEF': 0.45,  # (dim) daylight air temperature coefficient

          # parameters for the snowpack algorithm
          'SNOW_TCRIT': -6.0,   # (deg C) critical temperature for snowmelt
          'SNOW_TRATE': 0.042,  # (cm/degC/day) snowmelt rate

          # parameters for the radiation algorithm
          # (dim) stands for dimensionless values
          'TBASE': 0.870,       # (dim) max inst. trans., 0m, nadir, dry atm
          'ABASE': -6.1e-5,     # (1/Pa) vapor pressure effect on
                                # transmittance
          'C': 1.5,             # (dim) radiation parameter
          'B0': 0.031,          # (dim) radiation parameter
          'B1': 0.201,          # (dim) radiation parameter
          'B2': 0.185,          # (dim) radiation parameter
          'RAIN_SCALAR': 0.75,  # (dim) correction to trans. for rain day
          'DIF_ALB': 0.6,       # (dim) diffuse albedo for horizon correction
          'SC_INT': 1.32,       # (MJ/m2/day) snow correction intercept
          'SC_SLOPE': 0.096,    # (MJ/m2/day/cm) snow correction slope
          'site_elev': 0.,
          'base_elev': 0.,
          't_max_lr': 0.0065,
          't_min_lr': 0.0065,
          'site_isoh': None,
          'base_isoh': None,
          'site_lat': 0.,
          'site_slope': 0.,
          'site_aspect': 0.,
          'site_east_horiz': 0.,
          'site_west_horiz': 0.
         }

OPTIONS = {
           'SW_PREC_THRESH': 0.,
           'VP_ITER': 'VP_ITER_ALWAYS',
           'MTCLIM_SWE_CORR': False,
           'LW_CLOUD': 'CLOUD_DEARDORFF',
           'LW_TYPE': 'PRATA'
          }


def update(io_params: dict):
    """ 
    Update the configuration based on what was read in
    from the input file or config file.
    """
    # TODO: This should take into account if there are 
    #       keys named 'options' and 'consts' and put
    #       them in the right things.
    PARAMS.update(io_params)


from metsim.methods import mtclim

METHODS = {
           'mtclim': mtclim
          }



