"""
Stores default constants
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

try:
    # prefer to just use the dictionary stored in netCDF
    from netCDF4 import default_fillvals as FILL_VALUES
except ImportError:
    # but we can use these if we can't import netCDF4
    FILL_VALUES = {'S1': '\x00',
                   'f4': 9.969209968386869e+36,
                   'f8': 9.969209968386869e+36,
                   'i1': -127,
                   'i2': -32767,
                   'i4': -2147483647,
                   'i8': -9223372036854775806,
                   'u1': 255,
                   'u2': 65535,
                   'u4': 4294967295,
                   'u8': 18446744073709551614}

DEG_PER_REV = 360.0       # Number of degrees in full revolution
SEC_PER_RAD = 13750.9871  # seconds per radian of hour angle
RAD_PER_DAY = 0.017214    # radians of Earth orbit per julian day
RAD_PER_DEG = 0.01745329  # radians per degree
MIN_DECL = -0.4092797     # minimum declination (radians)
DAYS_OFF = 11.25          # julian day offset of winter solstice
# Note:  Make sure that 3600 % SW_RAD_DT == 0
SW_RAD_DT = 30.0          # timestep for radiation routine (seconds)

MA = 28.9644e-3         # (kg mol-1) molecular weight of air
R = 8.3143              # (m3 Pa mol-1 K-1) gas law constant
R_DRY = 287             # (J / degC * kg) Gas constant of dry air
G_STD = 9.80665         # (m s-2) standard gravitational accel.
P_STD = 101325.0        # (Pa) standard pressure at 0. m elevation
T_STD = 288.15          # (K) standard temp at 0. m elevation
CP = 1010.0             # (J kg-1 K-1) specific heat of air

KELVIN = 273.15
EPS = 0.62196351
DAYS_PER_YEAR = 365.25
SEC_PER_DAY = 86400
SEC_PER_MIN = 60
MIN_PER_HOUR = 60
HOURS_PER_DAY = 24
MIN_PER_DAY = 1440
SEC_PER_HOUR = 3600
STEFAN_B = 5.669e-8      # (W m^-22 K^-4) Stefan Boltzmann constant

M_PER_KM = 1000.0
MBAR_PER_BAR = 1000.0
MM_PER_CM = 10.0
MAX_PERCENT = 100.0

# parameters for the radiation algorithm
# (dim) stands for dimensionless values
C = 1.5             # (dim) radiation parameter
B0 = 0.031          # (dim) radiation parameter
B1 = 0.201          # (dim) radiation parameter
B2 = 0.185          # (dim) radiation parameter
TBASE = 0.870       # (dim) max inst. trans. 0m nadir dry atm
ABASE = -6.1e-5     # (1/Pa) vapor pressure effect on transmittance
SC_INT = 1.32       # (MJ/m2/day) snow correction intercept
SC_SLOPE = 0.096    # (MJ/m2/day/cm) snow correction slope
