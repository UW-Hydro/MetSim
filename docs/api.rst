.. currentmodule:: metsim

#############
API reference
#############

This page provides an auto-generated summary of metsim's API. For more details
and examples, refer to the relevant chapters in the main part of the
documentation.

MetSim
=======

Creating a MetSim Object
------------------------

.. autosummary::
   :toctree: generated/

   MetSim

Attributes
----------

.. autosummary::
   :toctree: generated/

   metsim.MetSim.methods
   metsim.MetSim.params
   metsim.MetSim.output
   metsim.MetSim.met_data
   metsim.MetSim.domain

   metsim.MetSim.load
   metsim.MetSim.launch
   metsim.MetSim.run
   metsim.MetSim.update
   metsim.MetSim.write
   metsim.MetSim.read

Methods
=======

Physics
-------

.. autosummary::
   :toctree: generated/

    metsim.physics.calc_pet
    metsim.physics.atm_pres
    metsim.physics.svp
    metsim.physics.svp_slope
    metsim.physics.solar_geom

MtClim
------

.. autosummary::
   :toctree: generated/

    metsim.methods.mtclim.run
    metsim.methods.mtclim.calc_t_air
    metsim.methods.mtclim.calc_prec
    metsim.methods.mtclim.calc_snowpack
    metsim.methods.mtclim.calc_srad_hum
    metsim.methods.mtclim.sw_hum_iter

Disagg
------
.. autosummary::
   :toctree: generated/

    metsim.disaggregate.disaggregate
    metsim.disaggregate.set_min_max_hour
    metsim.disaggregate.temp
    metsim.disaggregate.prec
    metsim.disaggregate.wind
    metsim.disaggregate.relative_humidity
    metsim.disaggregate.vapor_pressure
    metsim.disaggregate.longwave
    metsim.disaggregate.shortwave
