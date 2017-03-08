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
(number of minutes in a day), no disaggregation will occur. This value should be
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
exist it will be created.

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
algorithm. Default to ``0``.

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

domain_vars section
-------------------

