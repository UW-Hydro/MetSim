.. _configuration:

Configuration Specifications
============================
This page documents the various options and
parameters that can be set in the configuration
file.

MetSim Section
--------------

**Required Variables**

``time_step :: int``: The timestep to disaggregate in minutes.  If given as 1440
(number of minutes in a day), no disaggregation will occur. This value must
divide 1440 evenly.

``start :: str``: The time to start simulation given in the format 
``yyyy/mm/dd`` or ``yyyy/mm/dd:hh``.

``stop :: str``: The time to end simulation given in the format
``yyyy/mm/dd``.

``forcing :: path``: The path to the input forcing file(s).  See the section 
on __in_vars__ for more details.

``domain :: path``: The path to the input domain file.  See the section on 
__domain_vars__ for more details.

``out_dir :: path``: The location to write output to.  If this path doesn't 
exist, it will be created.

``in_format :: str``: A string representing the type of input files specified in
the ``forcing`` entry.  Can be one of the following: ``ascii``, ``binary``, 
``netcdf``.

``out_format :: str``: A string representing the type of output to write to 
``out_dir``.  Can be either ``netcdf`` or ``ascii``.

``method :: str``: A string representing the simulation methods to use.  The
current implementation only supports ``mtclim``.

``t_max_lr :: float``: The lapse rate (K/m) to be applied to the maximum daily 
temperature.

``t_min_lr :: float``: The lapse rate (K/m) to be applied to the minimum daily 
temperature.

**Optional Variables**

``out_prefix :: str``: The output file base name. Defaults to ``forcing``.

``verbose :: bool``: Whether to print output to ``stdout``.  Should be set using
the ``-v`` flag for command line usage.  This can be set for scripting purposes,
if desired. Set to ``1`` to print output; defaults to ``0``.

``base_elev :: float``: Base offset of elevation (m).  Defaults to ``0``.

``site_isoh :: float``: Isohyetal level offset for adjusting precipitation 
inputs.  Default to ``1``.

``base_isoh :: float``: Isoyetal base level for adjusting precipitation inputs. 
Defaults to ``1``.

``sw_prec_thresh :: float``: Minimum precipitation threshold to take into 
account when simulating incoming shortwave radiation.  Defaults to ``0``.

``mtclim_swe_corr :: bool``: Whether to activate MTCLIM's SWE correction 
algorithm. Default to ``False``.

``lw_cloud :: str``: Type of cloud correction to longwave radiation to apply. 
Can be either ``DEFAULT`` or ``CLOUD_DEARDORFF``.  Defaults to 
``CLOUD_DEARDORFF``.  Capitalization does not matter.

``lw_type :: str``: Type of longwave radiation parameterization to apply. Can be
one of the following: ``DEFAULT``, ``TVA``, ``ANDERSON``, ``BRUTSAERT``, 
``SATTERLUND``, ``IDSO``, or ``PRATA``.  Defaults to ``PRATA``.  Capitalization 
does not matter.

``tdew_tol :: float``: Convergence criteria for the iterative calculation of 
dewpoint temperature in MTCLIM.  Defaults to ``1e-3``.  

``tmax_daylength_fraction :: float`` : Weight for calculation of time of maximum
daily temperature.  Must be between ``0`` and ``1``.  Defaults to ``0.67``.

``out_vars :: list`` : List of variables to write to output.  Should be a list 
containing valid variables.  The list of valid variables is dependent on which 
simulation method is used, as well as whether disaggregation is used. Defaults 
to ``['temp', 'prec', 'shortwave', 'longwave', 'vapor_pressure', 'red_humid']``.


in_vars section
---------------
The ``in_vars`` section is where you can specify which variables are in your 
input data, and the corresponding symbols which MetSim will recognize. The 
format of this section depends on the value given in the ``in_format`` entry in 
the ``MetSim`` section of the configuration file.  See below for conventions for
each input type.


netcdf
``````
The ``in_vars`` section for NetCDF input acts as a mapping between the variable 
names in the input dataset to the variable names expected by MetSim.  The format
is given as ``netcdf_varname = metsim_varname``.  The minimum required variables 
given have ``metsim_varname``s corresponding to ``t_min``, ``t_max``, and 
``prec``; these variable names correspond to minimum daily temperature (Celcius), 
maximum daily temperature (Celcius), and precipitation (mm).

ascii
`````
The ``in_vars`` section for ASCII input acts similarly to the NetCDF input 
format, except for one key point.  Variables should be given as a tautology: the 
format is given as ``metsim_varname = metsim_varname``.  The order that the 
variables are given corresponds to the column numbers that they appear in the 
input files.  The minimum required variables are ``t_min``, ``t_max``, and 
``prec``; these variable names correspond to minimum daily temperature (Celcius), 
maximum daily temperature (Celcius), and precipitation (mm).

binary
``````
This section has an input style for binary files that mimics the VIC version 4 
input style.  Each line is specified as ``varname = scale cdatatype``, where 
``varname`` is the name that MetSim should use for the column, ``scale`` is a 
floating point scaling factor that should be applied after conversion from 
binary to floating point; the conversion applied by the ``scale`` is applied 
after the value in the input is converted from binary to the ``cdatatype`` 
specified for each variable.  Valid ``cdatatype``s are ``signed`` and 
``unsigned``.  ``signed`` values are interpreted as values which can be positive
or negative, whereas ``unsigned`` values are interpreted as values that can only
be greater than or equal to zero.

domain_vars section
-------------------
The ``domain_vars`` section is where information about the domain file is given.
Since the domain file is given as a NetCDF file this section has a similar 
format to that of the NetCDF input file format described above.  That is, 
entries should be of the form ``netcdf_varname = metsim_varname``. The minimum 
required variables have ``metsim_varname``s corresponding to ``lat``, ``lon``,
``mask``, and ``elev``; these variable names correspond to latitude, longitude, 
a mask of valid cells in the domain, and the elevation given in meters.
