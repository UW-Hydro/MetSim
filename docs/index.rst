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

.. _MtClim: http://www.ntsg.umt.edu/project/mt-clim.php
.. _VIC: https://github.com/UW-Hydro/VIC

MetSim consists of 3 main modules that govern the operation of 3
major aspects of its operation:

**1. Management of dataset preprocessing and IO**

The MetSim object provides high level support for setting up jobs
and infrastructure for running simulation/disaggregation
steps. It is the main interface through which the other modules
are accessed.

**2. Simulation of daily meteorological forcings**

The base implementation of the meteorological simulator is
based off of the algorithms described in [1]_. This component
has been designed to be flexible in allowing for alternative
implementations which may be specified during the setup of the
MetSim object.  The default implementation allows for the
daily simulation of:

 * Mean daily temperature
 * Incoming shortwave radiation
 * Cloud cover fraction
 * Potential evapotranspiration
 * Vapor pressure

**3. Disaggregation of daily simulation values to sub-daily timesteps**

Daily data from given input or simulated via the forcings generation
component of MetSim can be disaggregated down to sub-daily values at
intervals specified in minutes (provided they divide evenly into 24
hours).  The operation of these algorithms is also described in [1]_.
The variables estimated are:

 * Temperature
 * Vapor pressure
 * Relative and specific humidity
 * Air pressure
 * Cloud cover fraction
 * Longwave radiation
 * Shortwave radiation
 * Precipitation
 * Wind speed

For the "triangle" and "mix" methods of precipitation disaggregation,
doumentation can be found `here <PtriangleMethod.pdf>`_. This will eventually
be superceded by a journal article that is currently in review [7]_.

If you don't find what you're looking for here, check out MetSim's Github page.

Getting Started
===============
A tutorial for running MetSim and working with input/output data can be run
via binder here: https://github.com/UW-Hydro/MetSim-tutorial

Installation
------------
MetSim itself is a pure Python package, but its dependencies are not. You should
ensure that you have all of the required dependencies:

- Python 3.5 or greater
- `xarray <http://xarray.pydata.org/>`__ (0.10.9 or later)
- `pandas <http://pandas.pydata.org/>`__ (0.19.0 or later)
- `numba <http://numba.pydata.org/>`__ (0.31.0 or later)
- `netCDF4 <https://github.com/Unidata/netcdf4-python>`__
- `scipy <http://scipy.org/>`__


Then, install MetSim with pip or conda::

    $ pip install metsim

or::

    $ conda install -c conda-forge metsim

Alternatively, you can install MetSim directly from the source if you desire to::

    $ git clone https://github.com/UW-Hydro/MetSim.git
    $ cd MetSim
    $ python setup.py install
    $ py.test --verbose


Basic Usage
-----------
MetSim provides a simple command line interface which is primarily operated via
configuration files.  For more information about the options available to be set
in the configuration files see the `configuration <configuration.rst>`_ page.

Once installed, MetSim can be used from the command line via:

.. code-block:: bash

    usage: ms [-h] [-n NUM_WORKERS] [-s SCHEDULER] [-v] [--version] config

    positional arguments:
    config                Input configuration file

    optional arguments:
    -h, --help            show this help message and exit
    -n NUM_WORKERS, --num_workers NUM_WORKERS
                            Parallel mode: number of processes to use
    -s SCHEDULER, --scheduler SCHEDULER
                            Dask scheduler to use
    -v, --verbose         Increase the verbosity of MetSim
    --version             Name and version number

Bracketed flags are optional; ``-v`` activates verbose mode to print messages
about the status of a run, and ``-n`` activates parallelism.  The number given
after the ``-n`` flag is the number of processes to run. A good rule of thumb is
to use one less process than the number of processsors (or threads) that the
machine you are running on has.

.. warning::
    Users in environments where OpenMP is available may experience
    over-utilization of CPU resources, leading to lower performance. If you experience
    this issue try setting the `OMP_NUM_THREADS` environment variable to 1 before running
    MetSim.. This can be done in bash and similar shells by running
    `export OMP_NUM_THREADS=1`.


References
==========

.. [1] Bohn, T. J., B. Livneh, J. W. Oyler, S. W. Running, B. Nijssen, and D. P.
       Lettenmaier, 2013. Global evaluation of MTCLIM and related algorithms for
       forcing of ecological and hydrological models, Agricultural and Forest
       Meteorology, 176:38-49, doi:10.1016/j.agrformet.2013.03.003.

.. [2] Bristow, K.L., and G.S. Campbell, 1984. On the relationship between
       incoming solar radiation and daily maximum and minimum temperature.
       Agricultural and Forest Meteorology, 31:159-166.

.. [3] Running, S.W., R.R. Nemani, and R.D. Hungerford, 1987. Extrapolation of
       synoptic meteorological data in mountainous terrain and its use for
       simulating forest evaporation and photosynthesis. Canadian Journal of
       Forest Research, 17:472-483.

.. [4] Glassy, J.M., and S.W. Running, 1994. Validating diurnal climatology of
       the MT-CLIM model across a climatic gradient in Oregon. Ecological
       Applications, 4(2):248-257.

.. [5] Kimball, J.S., S.W. Running, and R. Nemani, 1997. An improved method for
       estimating surface humidity from daily minimum temperature. Agricultural
       and Forest Meteorology, 85:87-98.

.. [6] Thornton, P.E., and S.W. Running, 1999. An improved algorithm for
       estimating incident daily solar radiation from measurements of
       temperature, humidity, and precipitation. Agricultural and Forest
       Meteorology, 93:211-228.

.. [7] Bohn, T. J., K. M. Whitney, G. Mascaro, and E. R. Vivoni, 2019. A
       deterministic approach for approximating the diurnal cycle of
       precipitation for large-scale hydrological simulations. Journal of
       Hydrometeorology, 20(2):297-317. doi: 10.1175/JHM-D-18-0203.1.

Sitemap
=======
.. toctree::
    :maxdepth: 3

    examples
    data
    configuration
    configuration_ini
    api
    whats-new
