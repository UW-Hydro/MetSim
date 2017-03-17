.. metsim documentation master file, created by
   sphinx-quickstart on Mon Feb 27 17:34:14 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. _index:



METSIM: Meteorology Simulator
==================================
MetSim is a meteorological simulator and forcing disaggregator
for hydrologic modeling and climate applications.
Metsim is based on MtClim_ and the preprocessor from version 4 of the
VIC_ hydrologic model.

.. _MtClim: http://www.ntsg.umt.edu/project/mtclim
.. _VIC: https://github.com/UW-Hydro/VIC

MetSim consists of 3 main modules that govern the operation of 3 
major aspects of its operation:

**1. Management of dataset preprocessing and IO**

The MetSim object provides high level support for setting up jobs
and infrastructure for running simulation/disaggregation
steps. It is the main interface through which the other modules
are accessed. 

**2. Simulation of meteorological forcings**

The base implementation of the meteorological simulator is
based off of the algorithms described in [1]_. This component
has been designed to be flexible in allowing for alternative 
implementations which may be specified during the setup of the
MetSim object.  The default implementation allows for the 
daily simulation of:

 * Mean daily temperature
 * Snow water equivalent (SWE)
 * Incoming shortwave radiation
 * Cloud cover fraction
 * Potential evapotranspiration
 * Vapor pressure

**3. Disaggregation of daily simulation values to sub-daily timesteps**

Daily data from given input or simulated via the forcings generation
component of MetSim can be disaggregated down to sub-daily values at
intervals specified in minutes (provided they divide evenly into 24 
hours).  The operation of these algorithms is also described in [1]_.


Depending on the setup of these various components MetSim can be used
for different purposes.  For examples on how you can use MetSim 
see the :ref:`examples` page.

This documentation is a work in progress.
If you don't find what you're looking for here, check out MetSim's Github page.  

Getting Started
===============

Installation
------------
MetSim itself is a pure Python package, but its dependencies are not. The
easiest way to get everything installed is to use conda_. To install MetSim 
with its recommended dependencies using the ``conda`` command line tool::

    $ conda install metsim

.. _conda: http://conda.io/

We recommend using the community maintained `conda-forge 
<https://conda-forge.github.io/>`__ channel if you need difficult\-to\-build 
dependencies such as ``numba``::

    $ conda install -c conda-forge metsim

New releases may also appear in conda-forge before being updated in the default
channel.

If you are not using ``conda`` to install MetSim you should ensure that you have
all of the required dependencies:

- Python 3.5 or 3.6
- `xarray <http://xarray.pydata.org/>`__ (0.9.1 or later)
- `pandas <http://pandas.pydata.org/>`__ (0.19.0 or later)
- `numba <http://numba.pydata.org/>`__ (0.31.0 or later)
- `netCDF4 <https://github.com/Unidata/netcdf4-python>`__ 
- `scipy <http://scipy.org/>`__


If you don't use conda, be sure you have the required dependencies (numba and 
xarray) installed first. Then, install MetSim with pip::

    $ pip install metsim

To run the test suite after installing MetSim, install
`py.test <https://pytest.org>`__ (``pip install pytest``) and run
``py.test --verbose``.

Finally, you can install MetSim directly from the source if you desire to::

    $ git clone https://github.com/UW-Hydro/MetSim.git
    $ cd MetSim
    $ python setup.py install
    $ py.test --verbose


Basic Usage
-----------
MetSim provides a simple command line interface which is primarily operated via 
configuration files.  For more information about the options available to be set
in the configuration files see the :ref:`configuration` page.

Once installed, MetSim can be used from the command line via:

``ms /path/to/configuration [-v] [-n #]``

Bracketed flags are optional; ``-v`` activates verbose mode to print messages 
about the status of a run, and ``-n`` activates parallelism.  The number given 
after the ``-n`` flag is the number of processes to run. A good rule of thumb is
to use one less process than the number of processsors (or threads) that the 
machine you are running on has.

References
==========

.. [1] Bohn, T. J., B. Livneh, J. W. Oyler, S. W. Running, B. Nijssen, and D. P.
    Lettenmaier, 2013a: Global evaluation of MTCLIM and related algorithms for 
    forcing of ecological and hydrological models, Agr. Forest. Meteorol., 176, 
    38-49, doi:10.1016/j.agrformet.2013.03.003.

Sitemap
=======
.. toctree::
    :maxdepth: 3
 
    index
    examples
    configuration
    api
