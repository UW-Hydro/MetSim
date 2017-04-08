.. _data:

Output Specifications
=====================
The output variables that are available are dependent on the time step being used.  There are two cases:

Daily Output
------------

When ``time_step`` is set to 1440 in the configuration file, daily values are generated.  The following variables are available for output at a daily time step:

* ``t_min`` : Minimum temperature (also a required input value) (C)
* ``t_max`` : Maximum temperature (also a required input value) (C)
* ``prec`` : Precipitation (also a required input value) (mm/day)
* ``swe`` : Snow water equivalent (mm)
* ``vapor_pressure`` : Vapor pressure (kPa)
* ``swrad`` : Shortwave radiation (W/m^2)
* ``tskc`` : Cloud cover fraction
* ``pet`` : Potential evapotranpiration (mm/day)

Sub-daily Output
----------------

* ``temperature`` : Temperature (C)
* ``shortwave`` : Shortwave radiation (W/m^2)
* ``vapor_pressure`` : Vapor pressure (kPa)
* ``rel_humid`` : Relative humidity
* ``longwave`` : Longwave radiation (W/m^2)
* ``tsck`` : Cloud cover fraction
* ``wind`` : Wind speed (only if given as an input) (m/s)
