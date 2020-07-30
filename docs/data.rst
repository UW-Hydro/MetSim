.. _data:

Input Specifications
====================
There are 3 required input data sources that must be specified in the
configuration file or dictionary. Note that it is possible for a single file to
be specified for all three sources, provided that it has all of the required
data. For examples of the data see the ``tests/data`` directory within the
MetSim code.

Input forcings
--------------
Specified as ``forcing`` in the configuration file. This can either be the path
to a NetCDF file, or the path to a directory containing ASCII or binary data (in
the VIC4 format). The input forcing data is used to provide the base forcing
variables. The required variable data is minimum daily temperature, maximum
daily temperature, and daily precipitation.

The variable names can be mapped via the configuration file in the ``forcing_vars``
section. For more information about how to set up your configuration file see
the `configuration <configuration.rst>`_ page.

Domain file
-----------
Specified as ``domain`` in the configuration file. The domain file provides
information about the domain MetSim is to be run over. It is required to be a
NetCDF file. The domain requires the following variables to be valid:

1. ``mask``: This provides information about which grid cells are valid to run
MetSim on. Values that specify grid cells which should be processed are
specified via a positive, finite number (one or greater). Cells which MetSim
should ignore can be given as 0 or ``NaN``.

It is important to ensure that all valid locations in ``mask`` have data in
``elev`` and any other variables.  Failure to ensure this will result in
errors during runtime.

2. ``elev``: This provides elevation data (in m) used for calculation of solar
geometry. It only needs to be given at sites which are marked to be processed
via the ``mask`` variable.

The next two variables are only needed if ``prec_type`` = ``triangle`` or
``mix`` in the input file:

3. ``dur``: This provides the climatological monthly storm event duration (in
minutes) used for disaggregating daily precipitation according to the
"triangle" method. Requires one value for each month (12).

4. ``t_pk``: This provides the climatological monthly time to storm peak (in
minutes starting from midnight) used for disaggregating daily precipitation to
sub-daily time scales using the "triangle" method. Requires one value for
each month (12).

For more information about the "triangle" method see
`this description <PtriangleMethod.pdf>`_. If you use this feature, please
cite Bohn et al. (2019) as listed in the `references <index.rst#id10>`_.

A domain file for the CONUS+Mexico domain, at 0.0625 degree resolution, and
containing ``dur`` and ``t_pk`` values, is available `here
<https://zenodo.org/record/1402223#.XEZC4M2IZPY>`_.

State file
----------
The state file provides information about the history of each of the grid cells
to be processed. There are four required variables.

The first two are daily minimum and daily maximum temperatures for the 90 days
preceeding the start date specified in the configuration file.  They should be
specified as ``t_min`` and ``t_max`` respectively. Similarly precipitation
should be given as ``prec``.  These variables are used to generate seasonal
averages which are used in the calculation of shortwave and longwave radiation.

Output Specifications
=====================
.. ATTENTION::
    The ``time`` coordinate in MetSim's output is local to the location of each cell unless the ``utc_offset`` is set to
    ``True``! This means that for a single time slice in the NetCDF file all locations along a parallel (same latitude)
    will have the same solar geometry at that time.

The output variables that are available are dependent on the time step being used.  There are two cases:

Daily Output
------------

When ``time_step`` is set to 1440 in the configuration file, daily values are
generated. The following variables are available for output at a daily time
step:

* ``t_min`` : Minimum temperature (also a required input value) (C)
* ``t_max`` : Maximum temperature (also a required input value) (C)
* ``prec`` : Precipitation (also a required input value) (mm/day)
* ``vapor_pressure`` : Vapor pressure (Pa)
* ``shortwave`` : Shortwave radiation (W/m^2)
* ``tskc`` : Cloud cover fraction
* ``pet`` : Potential evapotranpiration (mm/day)
* ``wind`` : Wind speed (only if given as an input) (m/s)

Sub-daily Output
----------------
If disaggregation to shorter time steps is configured, the data is passed from the meteorology simulation module to the disaggregation module.
MetSim implements several variable-specific disaggregation routines.
Here we briefly mention the disaggregation procedures for completeness, but no substantial changes were made to the earlier algorithms.

Shortwave is disaggregated by multiplying the total daily shortwave by the fraction of radiation received in a given timestep (provided by the solar geometry module).
This calculation is corrected for cloud cover by assuming constant transmissivity throughout the day (which is calculated in the meteorological simulation module).
Temperature is disaggregated by estimating the time at which the daily maximum and daily minimum temperatures occur.
These are chosen so that the daily minimum temperature occurs at sunrise and the daily maximum temperature occurs at a fixed time during the day (which is configurable by the user as a parameter in the configuration file if desired).
Then a Hermite polynomial interpolation is used to obtain the full temperature timeseries at sub-daily time steps.
Vapor pressure is disaggregated by linearly interpolating between the saturation vapor pressure values calculated based on the daily minimum temperature and that are assumed to occur at the time of the daily minimum temperature.
An additional correction is made to ensure that the vapor pressure at any given time step does not exceed the saturation vapor pressure, which is calculated directly from the disaggregated temperature timeseries.
Air pressure is disaggregated by using the disaggregated temperature as well as the elevation data provided by the domain file.
Both specific and relative humidity are then disaggregated using the disaggregated temperature and air pressure time series.
If provided, wind speed is disaggregated, but is assumed to be constant throughout the day.

As part of the model configuration, the user can select from a number of different algorithms to estimate longwave radiation.
Sub-daily values are calculated with the selected method using the disaggregated values for vapor pressure and temperature.

Precipitation can be disaggregated in one of two ways.
The first and simplest way is to evenly spread the daily precipitation across the sub-daily time steps.
The second method requires two additional parameters to be specified in the domain file to represent the average precipitation duration and the time of peak precipitation for each cell.
The method then disaggregates precipitation by constructing a triangular kernel with total area equal to the daily precipitation centered at the time which is specified as the time of peak precipitation.

The following variables (and name conventions) are allowed to be specified for subdaily output:

* ``temp`` : Temperature (C)
* ``prec`` : Precipitation (mm/timestep)
* ``shortwave`` : Shortwave radiation (W/m^2)
* ``vapor_pressure`` : Vapor pressure (Pa)
* ``air_pressure`` : Air pressure (kPa)
* ``rel_humid`` : Relative humidity
* ``spec_humid`` : Specific humidity
* ``longwave`` : Longwave radiation (W/m^2)
* ``tsck`` : Cloud cover fraction
* ``wind`` : Wind speed (only if given as an input) (m/s)
