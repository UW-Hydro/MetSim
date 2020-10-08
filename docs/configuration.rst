.. _configuration:

Configuration Specifications
============================
This page documents the various options and
parameters that can be set in the configuration
file. An example configuration file can be found on the examples page.

MetSim Section
--------------

**Required Variables**

``time_step :: int``: The timestep to disaggregate in minutes.  If given as 1440
(number of minutes in a day), no disaggregation will occur. This value must
divide 1440 evenly.

``start :: str``: The time to start simulation given in the format
``yyyy/mm/dd``

``stop :: str``: The time to end simulation given in the format
``yyyy/mm/dd``.

``forcing :: path``: The path to the input forcing file(s).  See the section
on __forcing_vars__ for more details.

``domain :: path``: The path to the input domain file.  See the section on
__domain_vars__ for more details.

``state :: path``: The path to the input state file.

``out_dir :: path``: The location to write output to.  If this path doesn't
exist, it will be created.

``forcing_fmt :: str``: A string representing the type of input files specified in
the ``forcing`` entry.  Can be one of the following: ``ascii``, ``binary``,
``netcdf``, or ``data``.

**Optional Variables**

``method ::str``: The method to use for estimation of meteorological quantities.
This can be either ``mtclim`` to estimate missing variables or ``passthrough`` if
some of the meteorological variables have already been estimated (for example, by
DayMet, PRISM, or GridMET). Defaults to ``mtclim``.

``out_prefix :: str``: The output file base name. Defaults to ``forcing``.

``out_precision :: str``: Precision to use when writing output.  Defaults to
``f8``.  Can be either ``f4`` or ``f8``.

``time_grouper :: str``: Whether to chunk up the timeseries into pieces for
processing. This option is useful to set for when you are limited on
memory.  Each chunk of output is written as ``{out_prefix}_{date_range}`` when
active. Any valid ``pandas.TimeGrouper`` string may be used (e.g. use '10AS'
for 10 year chunks).

``verbose :: bool``: Whether to print output to ``stdout``.  Should be set using
the ``-v`` flag for command line usage.  This can be set for scripting purposes,
if desired. Set to ``1`` to print output; defaults to ``0``.

``sw_prec_thresh :: float``: Minimum precipitation threshold to take into
account when simulating incoming shortwave radiation.  Defaults to ``0``.

``rain_scalar :: float``: Scale factor for calculation of cloudy sky
transmittance.  Defaults to ``0.75``, range should be between ``0`` and
``1``.

``utc_offset :: bool``: Whether to use UTC timecode offsets for shifting
timeseries. Without this option all times should be considered local to
the gridcell being processed. Large domain runs probably want to set this
option to ``True``.

``lw_cloud :: str``: Type of cloud correction to longwave radiation to apply.
Can be either ``DEFAULT`` or ``CLOUD_DEARDORFF``.  Defaults to
``CLOUD_DEARDORFF``.  Capitalization does not matter.

``lw_type :: str``: Type of longwave radiation parameterization to apply. Can be
one of the following: ``DEFAULT``, ``TVA``, ``ANDERSON``, ``BRUTSAERT``,
``SATTERLUND``, ``IDSO``, or ``PRATA``.  Defaults to ``PRATA``.  Capitalization
does not matter.

``tdew_tol :: float``: Convergence criteria for the iterative calculation of
dewpoint temperature in MtClim.  Defaults to ``1e-6``.

``tmax_daylength_fraction :: float`` : Weight for calculation of time of maximum
daily temperature.  Must be between ``0`` and ``1``.  Defaults to ``0.67``.

``tday_coef :: float``: Scale factor for calculation of daily mean temperature.
Defaults to ``0.45``, range should be between ``0`` and ``1``.

``lapse_rate :: float``: Used to calculate atmospheric pressure. Defaults to
``0.0065 K/m``.

``out_vars :: list`` : List of variables to write to output.  Should be a list
containing valid variables.  The list of valid variables is dependent on which
simulation method is used, as well as whether disaggregation is used. Defaults
to ``['temp', 'prec', 'shortwave', 'longwave', 'vapor_pressure', 'red_humid']``.

``prec_type :: str``: Type of precipitation disaggregation method to use. Can be
one of the following: ``uniform``, ``triangle``, or ``mix``. Defaults to
``uniform``.  Capitalization does not matter. Under ``uniform`` method,
precipitation is disaggregated by dividing uniformly over all sub-daily
timesteps. Under ``triangle`` the "triangle" method is employed whereby daily
precipitation is distributed assuming an isosceles triangle shape with peak and
width determined from two domain variables, ``t_pk`` and ``dur``.  Under
``mix``, the "uniform" method is used on days when ``t_min`` < 0 C, and
"triangle" is used on all other days; this hybrid method retains the improved
accuracy of "triangle" in terms of warm season runoff but avoids the biases
in snow accumulation that the "triangle" method sometimes yields due to fixed
event timing within the diurnal cycle of temperature. A domain file for the
CONUS+Mexico domain, containing the ``dur`` and ``t_pk`` parameters is
available at: `<https://zenodo.org/record/1402223#.XEI-mM2IZPY>`.  For more
information about the "triangle" method see :doc:`PtriangleMethod.pdf`.

For more information about input and output variables see the :ref:`data` page.
::

    # Comments begin with hashtags
    # The first non-comment line must begin with the following:
    Metsim:
        time_step: int
        start: YYYY-MM-DD
        stop: YYYY-MM-DD

        # Paths to input files
        forcing: str
        domain: str
        state: str

        # Output file specification
        out_dir: str
        out_prefix: str

        # Algorithmic controls
        utc_offset: bool
        prec_type: str
        lw_type: str
        lw_cloud: str


chunks section
--------------
The ``chunks`` section describes how parallel computation should be grouped
in space. For example, to parallelize over 10 by 10 chunks of latitude and
longitude (with netcdf dimensions named ``lat`` and ``lon``, respectively) you would use:
::

    chunks:
        lat: 10
        lon: 10

Alternatively, for an HRU based run chunked into 50 element jobs you would use:
::

    chunks:
        hru: 50

As a general rule of thumb, try to evenly chunk the domain in such a way that
the number of jobs to run is some multiple of the number of processors you wish
to run on.

forcing_vars and state_vars section
---------------
The ``forcing_vars`` and ``state_vars`` sections are where you can specify which
variables are in your input data, and the corresponding symbols which MetSim will
recognize.  The ``in_vars`` section for acts as a mapping between the variable
names in the input dataset to the variable names expected by MetSim.  The format
is given as ``metsim_varname: netcdf_varname``.  The minimum required variables
given have ``metsim_varname``\s corresponding to ``t_min``, ``t_max``, and
``prec``; these variable names correspond to minimum daily temperature (Celcius),
maximum daily temperature (Celcius), and precipitation (mm/day).

domain_vars section
-------------------
The ``domain_vars`` section is where information about the domain file is given.
Since the domain file is given as a NetCDF file this section has a similar
format to that of the NetCDF input file format described above.  That is,
entries should be of the form ``metsim_varname = netcdfvarname``. The minimum
required variables have ``metsim_varname``\s corresponding to ``lat``, ``lon``,
``mask``, and ``elev``; these variable names correspond to latitude, longitude,
a mask of valid cells in the domain, and the elevation given in meters. If
``prec_type`` = ``triangle`` or ``mix``, two additonal variables are required
including ``dur`` and ``t_pk`` for disaggregating daily precipitation according
to the "triangle" method.

out_vars section
----------------
The ``out_vars`` section is where you can specify the output variables that you
want to include. There are two formats for this section. The first is the old format,
which we provide backwards compatibility for. You simply specify in the top level
``[MetSim]`` section a list of output variables with the names used by MetSim. They
will be written out with the same names used internally. Available options are
dependent on whether daily or subdaily output is being generated. Options for
daily output are:

- pet
- shortwave
- t_max
- t_min
- tskc

Options for subdaily output are:

 - prec
 - shortwave
 - longwave
 - temp
 - vapor_pressure
 - air_pressure
 - tskc
 - rel_humid
 - spec_humid
 - wind

The syntax for output specification is as follows:
::

    out_vars:
        metsim_varname:
            out_name: str
            units: str

unit conversions
================
The ``out_vars`` section allows for specification of some simple unit conversions
for MetSim output. The allowed options are as follows (invalid options will revert
to the default after issuing a warning):

 * prec
   - ``mm timestep-1`` (default)
   - ``mm s-1``
   - ``mm h-1``
 * pet (daily output only)
   - ``mm timestep-1`` (default)
   - ``mm s-1``
   - ``mm h-1``
 * t_max (daily output only)
   - ``C`` (default)
   - ``K``
 * t_min (daily output only)
   - ``C`` (default)
   - ``K``
 * temp
   - ``C`` (default)
   - ``K``
 * vapor_pressure
   - ``Pa`` (default)
   - ``hPa``
   - ``KPa``
 * air_pressure
   - ``kPa`` (default)
   - ``hPa``
   - ``Pa``
 * tskc (cloud fraction)
   - ``fraction`` (default)
   - ``%``
 * rel_humid
   - ``%`` (default)
   - ``fraction``

constant_vars section
-------------------
The ``constant_vars`` section is optional and allows you to set some of the
forcing inputs to a constant value. The specification simply consists of entries
of the form ``metsim_varname: value``, where ``value`` is a number that can be
converted to a double. There can only be one entry per line. If the
``metsim_varname`` corresponds to an entry that is already in the ``forcing_vars``
section, then the constant value will take precedence. In the current
implementation there must be at least one non-constant entry in ``forcings_vars``
(i.e. at least one entry that is not also in ``constant_vars``).

For example:
::

    constant_vars:
        wind: 2.0

will result in a constant wind field in the output file. In this case ``wind``
does not need to be specified in the ``forcing_vars`` section. If it was, it
will still be set to a constant value of 2 m/s.

Similarly:
::

    constant_vars:
        t_max = 30.0
        t_min = 10.0

will result in output with a diurnal cycle in which the temperature varies at
all locations between 10C and 30C. However, all estimation and disaggregation
routines are still evaluated, with constant ``t_max`` and ``t_min`` as input.
