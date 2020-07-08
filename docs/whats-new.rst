.. currentmodule:: MetSim

What's New
==========

.. _whats-new.2.3.0:

v2.3.0
------
Enchancements
~~~~~~~~~~~~~
- Allow for variable renaming in INI configuration files.
- Added capability for new YAML configuration file format
- Added capability simple unit conversions on output when
  using the YAML configuration file format

Bug fixes
~~~~~~~~~
- Fixed a bug where `utc_offset` causes radiation to be
  incorrectly scaled

.. _whats-new.2.2.0:

v2.2.2
-------
Bug fixes
~~~~~~~~~
- Fixed bug where `utc_offset` doesn't get converted to the
  correct boolean when reading the configuration.

v2.2.l
------
Bug fixes
~~~~~~~~~
- Fixed bug where timestamps got duplicated for sub-hourly
  disaggregation time periods
- Fixed bug where polar latitudes caused error in setting
  the rise and set times for temperature disaggregation

v2.2.0
-------
Enhancements
~~~~~~~~~~~~
- Can now specify ``period_ending`` in the configuration to move
  timestamps to end of period instead of beginning of period
- Addition of tutorial in README.md and main documentation
- Addition of paper to be submitted to JOSS

.. _whats-new.2.1.2:

v.2.1.2 (under development)
---------

Bug fixes
~~~~~~~~~
- Can now handle dimensions without coordinate variables, which
  previously caused a bug in the chunk selection of the worker
  processes. The fix is to simply add a sequential coordinate
  any time this occurs.

.. _whats-new.2.1.1:

v.2.1.1
-------
Enhancements
~~~~~~~~~~~~
- Allow for ``--version`` flag on command-line.

.. _whats-new.2.0.1:

v.2.0.1
-------
Enhancements
~~~~~~~~~~~~
- Allow for specification of constant fields, through addition
  of an optional ``constant_vars`` section.


.. _whats-new.2.0.0:

v.2.0.0
-------

Enhancements
~~~~~~~~~~~~
- Implemented UTC offsets, which puts all gridcell times in reference to UTC.
- Moved parallelism to dask, which allows for greater scalability and
  significantly less memory overhead.

Bug fixes
~~~~~~~~~
- Disallow timesteps > 6 hours, which raised errors.
- Raise error when t_min > t_max at beginning of runtime.

.. _whats-new.1.1.1:

v.1.1.1
-------

Enhancements
~~~~~~~~~~~~

- Added option to disaggregate precipitation via a triangular hyetograph.
  (:issue:`42`).
  By `Kristen Whitney <https://github.com/kwhitney727>`_ and `Theodore Bohn
  <https://github.com/tbohn>`_.

Bug fixes
~~~~~~~~~
- Fixed a bug where if `iter_dims` is not `[lat, lon]` the selected `lat` value
  that goes into `solar_geom` ends up as a list. The fix is also added for elevation
  and longitude, for redundancy.  Fixes :issue:`132`.
  By `Andrew Bennett <https://github.com/arbennett>`_.


.. _whats-new.1.1.0:

v1.1.0
------

Enhancements
~~~~~~~~~~~~

- Added option to use forcing start/stop dates to define run length (:issue:`93`).
  By `Joe Hamman <https://github.com/jhamman>`_.
- Added option a flexible time grouper when chunking MetSim runs (:issue:`93`).
  By `Joe Hamman <https://github.com/jhamman>`_.
- Improved configuration validation by checking for correctness of output variables (:issue:`96`)
  By `Andrew Bennett <https://github.com/arbennett>`_.
- Added option to skip reading ``swe`` variable from state file if it is not
  going to be used by MtClim. (:issue:`103`). By `Joe Hamman <https://github.com/jhamman>`_.
- Added support for supplying a glob-like file path or multiple input forcing
  files (netCDF) (:issue:`126`). By `Joe Hamman <https://github.com/jhamman>`_.
- Refactored ``mtclim`` and ``disaggregate`` functions to reduce interdependency and
  increase modularity. By `Andrew Bennett <https://github.com/arbennett>`_.
- Removed ``swe`` calculations. By `Andrew Bennett <https://github.com/arbennett>`_.

Bug fixes
~~~~~~~~~
- Fixed bug where output files were not written with the appropriate calendar
  encoding attribute (:issue:`97`).
  By `Joe Hamman <https://github.com/jhamman>`_.
- Fixed a bug where invalid timesteps were used in subdaily disaggregation.
  Added a clear error message explaining that subdaily timesteps must be evenly
  divisible into 24 hours and less than 6 hours in length. (:issue:`110`).
  By `Joe Hamman <https://github.com/jhamman>`_.
- Fixed a bug during disaggregation when ``t_min > t_max``.  This now raises
  an exception.
  By `Andrew Bennett <https://github.com/arbennett>`_.
