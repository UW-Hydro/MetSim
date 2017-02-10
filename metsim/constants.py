"""
Stores default constants, parameters, and options.
"""

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
G_STD = 9.80665         # (m s-2) standard gravitational accel.
P_STD = 101325.0        # (Pa) standard pressure at 0. m elevation
T_STD = 288.15          # (K) standard temp at 0. m elevation
CP = 1010.0             # (J kg-1 K-1) specific heat of air
LR_STD = 0.0065         # (-K m-1) standard temperature lapse rate

KELVIN = 273.15
EPS = 0.62196351
DAYS_PER_YEAR = 365.25
SEC_PER_DAY = 86400
MIN_PER_HOUR = 60
HOURS_PER_DAY = 24
SEC_PER_HOUR = 3600
STEFAN_B = 5.669e-8      # (W m^-22 K^-4) Stefan Boltzmann constant
TDAY_COEF = 0.45         # (dim) daylight air temperature coefficient

M_PER_KM = 1000.0
MBAR_PER_BAR = 1000.0
MM_PER_CM = 10.0
MAX_PERCENT = 100.0

# parameters for the snowpack algorithm (from mtclim)
SNOW_TCRIT = -6.0   # (deg C) critical temperature for snowmelt
SNOW_TRATE = 0.042  # (cm/degC/day) snowmelt rate

# parameters for the radiation algorithm
# (dim) stands for dimensionless values
TBASE = 0.870       # (dim) max inst. trans. 0m nadir dry atm
ABASE = -6.1e-5     # (1/Pa) vapor pressure effect on transmittance
C = 1.5             # (dim) radiation parameter
B0 = 0.031          # (dim) radiation parameter
B1 = 0.201          # (dim) radiation parameter
B2 = 0.185          # (dim) radiation parameter
RAIN_SCALAR = 0.75  # (dim) correction to trans. for rain day
DIF_ALB = 0.6       # (dim) diffuse albedo for horizon correction
SC_INT = 1.32       # (MJ/m2/day) snow correction intercept
SC_SLOPE = 0.096    # (MJ/m2/day/cm) snow correction slope
SW_PREC_THRESH = 0.
MTCLIM_SWE_CORR = False
VP_ITER = 'VP_ITER_ALWAYS'
LW_CLOUD = 'CLOUD_DEARDORFF'
LW_TYPE = 'PRATA'
TDEW_TOL = 1e-3
TMAX_DAYLENGTH_FRACTION = 0.67
