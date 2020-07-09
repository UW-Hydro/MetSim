.. _examples:

Examples
========

Basics
------
Provided in the source are several examples that can help you to
get started using MetSim. They are located in the ``examples``
directory.  For demonstration here is an example YAML configuration file::

    # This is an example of an input file for MetSim
    # Overall configuration, specification of parameters and input/output
    # paths goes in the "MetSim" section
    MetSim:
        # Time step in minutes
        time_step: 30
        # Forcings begin here (year-month-day)
        start: 1950-1-1
        # Forcings end at this date (year-month-day)
        stop: 1950-1-31
        # Input and output directories
        forcing: './metsim/data/test.nc'
        domain: './metsim/data/tiny_domain.nc'
        state: './metsim/data/state_nc.nc'
        forcing_fmt: 'netcdf'
        in_format: 'netcdf'
        out_dir: './results'
        out_prefix: 'yaml_output'
        prec_type: 'triangle'
        utc_offset: True

    out_vars:
        temp:
            out_name: 'airtemp'
            units: 'K'

        prec:
            out_name: 'pptrate'
            units: 'mm/s'

        shortwave:
            out_name: 'SWradAtm'

        spec_humid:
            out_name: 'spechum'

        air_pressure:
            out_name: 'airpres'
            units: 'kPa'

        wind:
            out_name: 'windspd'

    chunks:
        lat: 3
        lon: 3

    forcing_vars:
        # Format is metsim_name: input_name
        prec  : 'Prec'
        t_max : 'Tmax'
        t_min : 'Tmin'

    state_vars:
        # Format is metsim_name: input_name
        prec  : 'Prec'
        t_max : 'Tmax'
        t_min : 'Tmin'

    domain_vars:
        # Format is metsim_name: input_name
        lat  : 'lat'
        lon  : 'lon'
        mask : 'mask'
        elev : 'elev'
        t_pk : 't_pk'
        dur  : 'dur'

    constant_vars:
        wind : 2.0

This is a minimal configuration file for MetSim. For a complete description of the
input format see the `configuration <configuration.rst>`_ page.

To run this example from the command line, once you have installed
MetSim, use the following command:

``ms path/to/example.yaml --verbose``

This will run MetSim and disaggregate to half-hourly data, and write
out the results in a NetCDF file located in the directory specified
under ``out_dir`` in the configuration file (here ``./results``).
The addition of the ``--verbose`` flag provides some
information back to you as MetSim runs.  In the absence of this
flag MetSim will quietly run in the background until finished, or
some error has occurred.


Generating daily values
-----------------------
Daily values can be output by specifying a ``time_step`` of ``1440`` in the
configuration file, such as the one shown in the previous section. This will
prevent MetSim's disaggregation routines from being run, and the results written
out will be daily values.
