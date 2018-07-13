.. currentmodule:: MetSim

What's New
==========

.. _whats-new.1.1.1:

v.1.1.1
-------

Bug fixes
~~~~~~~~~
- Fixed a bug where if `iter_dims` is not `[lat, lon]` the selected `lat` value
  that goes into `solar_geom` ends up as a list. The fix is also added for elevation
  and longitude, for redundancy.  Fixes :issue:`132`.
  By `Andrew Bennett <https://github.com/arbennett>`


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
  By `Andrew Bennett <https://github.com/arbennett>`
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
