.. _examples:

Examples
========

Basics
------
Provided in the source are several examples that can help you to 
get started using using MetSim. They are located in the ``examples``
directory.  We will look at the ``example_nc.conf`` file.  It's 
contents are:

.. code-block:: ini
    
    # This is an example of an input file for MetSim
    [MetSim]
    
    # Time step in minutes
    time_step = 60
    
    # Forcings begin here (year/month/day:hour) (hour optional)
    start = 1950/1/1:0
    
    # Forcings end at this date (year/month/day)
    stop = 1950/1/31
    
    # Input and output directories
    forcing = ./tests/data/test.nc
    domain  = ./tests/data/domain.nc
    out_dir = ./results
    out_prefix = forcing
    in_format = netcdf
    out_format = netcdf
    
    # How to disaggregate
    method = mtclim
    
    t_max_lr = 0.0065
    t_min_lr = 0.0065
    
    
    [in_vars]
    Prec = prec
    Tmax = t_max
    Tmin = t_min
    
    [domain_vars]
    lat = lat
    lon = lon
    mask = mask
    elev = elev

This is a minimal configuration file for MetSim, and contains 
3 sections.  The first section, ``[MetSim]`` describes some
basic settings such as the locations of data and parameters
used in calculations.  For a complete description of the 
input format see :ref:`configuration`.  The key things to note
in this section are the ``forcing``, ``domain``, ``in_format``,
and ``out_format`` options.  The ``forcing`` and ``domain`` 
options refer to the two types of required input, and the ``in_format``
and ``out_format`` options tell MetSim how they should be treated.

The second section ``[in_vars]`` describes the variables in the
dataset provided in the ``forcing`` option of the first section.
The left side of the assignment is the name of the variable given
in the ``forcing`` dataset, while the right hand side is the 
name the variable should be given within MetSim.  Note that the
variables shown here are the minimum required set to run the
forcing generation. The names given on the right hand side are
also important to name correctly, as they are referenced internally.
If you are unsure what variable names are used internally see the 
:ref:`configuration` page for a full breakdown.

To run this example from the command line, once you have installed
MetSim, use the following command:

``ms path/to/example_nc.conf --verbose``

This will run MetSim and disaggregate to hourly data, and write 
out the results in a NetCDF file located in the directory specified
under ``out_dir`` in the configuration file (here ``./results``).
The addition of the ``--verbose`` flag should also provide some
information back to you as MetSim runs.  In the absence of this
flag MetSim will quietly run in the background until finished, or
some error has occurred.


Generating daily values
-----------------------

Generating hourly values
------------------------


Translating formats of daily values
-----------------------------------

.. warning:: This section `only` applies to daily input and output.


