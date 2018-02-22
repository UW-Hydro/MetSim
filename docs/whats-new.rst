.. currentmodule:: MetSim

What's New
==========

.. _whats-new.1.1.0:

v1.1.0 (unreleased)
-------------------

Enhancements
~~~~~~~~~~~~

- Added option to use forcing start/stop dates to define run length (:issue:`93`).
  By `Joe Hamman <https://github.com/jhamman>`_.
- Added option a flexible time grouper when chunking MetSim runs (:issue:`93`).
  By `Joe Hamman <https://github.com/jhamman>`_.
- Improved configuration validation by checking for correctness of output variables (:issue:`96`)
  By `Andrew Bennett <https://github.com/arbennett>`
- Added option to skip reading ``swe`` variable from state file if it is not
  going to be used by MtClim. (:issue:`XX`). By `Joe Hamman <https://github.com/jhamman>`_.

Bug fixes
~~~~~~~~~
- Fixed bug where output files were not written with the appropriate calendar
  encoding attribute (:issue:`97`).
  By `Joe Hamman <https://github.com/jhamman>`_.
- Fixed a bug where invalid timesteps were used in subdaily disaggregation.
  Added a clear error message explaining that subdaily timesteps must be evenly
  divisible into 24 hours and less than 6 hours in length. (:issue:`110`).
  By `Joe Hamman <https://github.com/jhamman>`_.
