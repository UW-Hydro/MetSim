.. metsim documentation master file, created by
   sphinx-quickstart on Mon Feb 27 17:34:14 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

METSIM: Meteorology Simulator
==================================

MetSim is a meteorological simulator and forcing disaggregator
for hydrologic modeling and climate applications.
Metsim is based on MtClim and the preprocessor from version 4 of the
VIC hydrologic model.

MetSim consists of 3 main modules that govern the operation of 3 
major aspects of its operation:

**Management of dataset preprocessing and IO**

The MetSim object provides high level support for setting up jobs
and infrastructure for running simulation/disaggregation
steps. It is the main interface through with the other modules
are accessed. 

**Simulation of meteorological forcings**

The base implementation of the meteorological simulator is
based off of the algorithms described in [1]_. This component
has been designed to be flexible in allowing for alternative 
implementations which may be specified during the setup of the
MetSim object.  The default implementation allows for the 
daily simulation of:

 * Mean daily temperature
 * SWE
 * Incoming shortwave radiation
 * Cloud cover fraction
 * Potential evapotranspiration
 * Vapor pressure

**Disaggregation of daily simulation values to sub-daily timesteps**

Daily data from given input or simulated via the forcings generation
component of MetSim can be disaggregaed down to sub-daily values at
intervals specified in minutes (provided they divide evenly into 24 
hours).  The operation of these algorithms is also described in [1]_.


Depending on the setup of these various components MetSim can be used
for different purposes.  For examples on how you can use MetSim 
see the :ref:`examples` page.

This documentation is a work in progress.
If you don't find what you're looking for here, check out metsim's Github page.  

References
==========

.. [1] Bohn, T. J., B. Livneh, J. W. Oyler, S. W. Running, B. Nijssen, and D. P.
    Lettenmaier, 2013a: Global evaluation of MTCLIM and related algorithms for 
    forcing of ecological and hydrological models, Agr. Forest. Meteorol., 176, 
    38-49, doi:10.1016/j.agrformet.2013.03.003.

.. toctree::
   :maxdepth: 1

   faq
   examples
   installing
   api
   whats-new

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
