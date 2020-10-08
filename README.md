METSIM: Meteorology Simulator
=============================
| MetSim Links & Badges              |                                                                             |
|------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| JOSS Paper             | [![DOI](https://joss.theoj.org/papers/10.21105/joss.02042/status.svg)](https://doi.org/10.21105/joss.02042) |
| MetSim Documentation      | [![Documentation Status](http://readthedocs.org/projects/metsim/badge/?version=develop)](http://metsim.readthedocs.io/en/develop/?badge=develop) |
| Travis-CI Build           | [![Build Status](https://travis-ci.org/UW-Hydro/MetSim.png)](https://travis-ci.org/UW-Hydro/MetSim) |
| License                | [![GitHub license](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://raw.githubusercontent.com/UW-Hydro/MetSim/master/LICENSE) |
| Current Release DOI    | [![DOI](https://zenodo.org/badge/69834400.svg)](https://zenodo.org/badge/latestdoi/69834400) |
| MetSim Tutorial        | https://github.com/UW-Hydro/MetSim-tutorial |

MetSim is a meteorological simulator and forcing disaggregator for
hydrologic modeling and climate applications. Metsim is based on
[MtClim](http://www.ntsg.umt.edu/project/mtclim)
and the preprocessor from version 4 of the [VIC hydrologic
model](https://github.com/UW-Hydro/VIC).

MetSim consists of 3 main modules that govern the operation of 3 major
aspects of its operation:

**1. Management of dataset preprocessing and IO**

The MetSim object provides high level support for setting up jobs and
infrastructure for running simulation/disaggregation steps. It is the
main interface through which the other modules are accessed.

**2. Simulation of daily meteorological forcings**

The base implementation of the meteorological simulator is based off of
the algorithms described in[1]. This component has been designed to be
flexible in allowing for alternative implementations which may be
specified during the setup of the MetSim object. The default
implementation allows for the daily simulation of:

-   Mean daily temperature
-   Incoming shortwave radiation
-   Cloud cover fraction
-   Potential evapotranspiration
-   Vapor pressure

For the "triangle" and "mix" methods of precipitation disaggregation,
doumentation can be found [here](https://github.com/UW-Hydro/MetSim/blob/develop/docs/PtriangleMethod.pdf).
This will eventually
be superceded by a journal article that is currently in review [7].

**3. Disaggregation of daily simulation values to sub-daily timesteps**

Daily data from given input or simulated via the forcings generation
component of MetSim can be disaggregated down to sub-daily values at
intervals specified in minutes (provided they divide evenly into 24
hours). The operation of these algorithms is also described in [1].
The variables estimated are:

-  Temperature
-  Vapor pressure
-  Relative and specific humidity
-  Air pressure
-  Cloud cover fraction
-  Longwave radiation
-  Shortwave radiation
-  Precipitation
-  Wind speed

Getting Started
===============
A tutorial for running MetSim and working with input/output data can be run
via binder here: https://github.com/UW-Hydro/MetSim-tutorial

Installation
------------

MetSim itself is a pure Python package, but its dependencies are not.
You should ensure that you have all of the required dependencies:

-   Python 3.5 or greater
-   [xarray](http://xarray.pydata.org/) (0.10.9 or later)
-   [pandas](http://pandas.pydata.org/) (0.19.0 or later)
-   [numba](http://numba.pydata.org/) (0.31.0 or later)
-   [netCDF4](https://github.com/Unidata/netcdf4-python)
-   [scipy](http://scipy.org/)

Then, install MetSim with pip or conda:

    $ pip install metsim

or

    $ conda install -c conda-forge metsim

Alternatively, you can install MetSim directly from the source if you desire to:

    $ git clone https://github.com/UW-Hydro/MetSim.git
    $ cd MetSim
    $ python setup.py install
    
If you are installing from source you may wish to also run the tests. 
You can do this from the MetSim directory with the command:

    $ pytest --verbose

Basic Usage
-----------

MetSim provides a simple command line interface which is primarily
operated via configuration files. For more information about the options
available to be set in the configuration files see the configuration
page in the full [documentation](http://metsim.readthedocs.io/en/develop/).


Once installed, MetSim can be used from the command line via:

`ms /path/to/configuration [-v] [-n #]`

Bracketed flags are optional; `-v` activates verbose mode to print
messages about the status of a run, and `-n` activates parallelism. The
number given after the `-n` flag is the number of processes to run. A
good rule of thumb is to use one less process than the number of
processsors (or threads) that the machine you are running on has.

:exclamation: Users in environments where OpenMP is available may experience
over-utilization of CPU resources, leading to lower performance. If you experience
this issue try setting the `OMP_NUM_THREADS` environment variable to 1 before running
MetSim.. This can be done in bash and similar shells by running
`export OMP_NUM_THREADS=1`.

Citing MetSim
=============
If you use MetSim in your work and would like to cite it you can use our JOSS paper as:

 > Bennett et al., (2020). MetSim: A Python package for estimation and disaggregation of meteorological data. Journal of Open Source Software, 5(47), 2042, https://doi.org/10.21105/joss.02042

Or in BibTeX:
``` 
 @article{Bennett2020,
  doi = {10.21105/joss.02042},
  url = {https://doi.org/10.21105/joss.02042},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {47},
  pages = {2042},
  author = {Andrew Bennett and Joseph Hamman and Bart Nijssen},
  title = {MetSim: A Python package for estimation and disaggregation of meteorological data},
  journal = {Journal of Open Source Software}
}
```

Acknowledgements
================
MetSim has greatly benefited from the user community, who have contributed code, tested features, provided feedback, and helped with documentation.
We would like to thank and acknowledge the work of Ted Bohn, Andy Wood, Kristen Whitney, Yifan Cheng, Liz Clark, Oriana Chegwidden, Ethan Gutmann, Kostas Andreadis, Thomas Remke, Ed Maurer, and Philipp Sommer for their help.


References
==========

[1]: Bohn, T. J., B. Livneh, J. W. Oyler, S. W. Running, B. Nijssen,
    and D. P. Lettenmaier, 2013a: Global evaluation of MTCLIM and
    related algorithms for forcing of ecological and hydrological
    models, Agr. Forest. Meteorol., 176, 38-49,
    <doi:10.1016/j.agrformet.2013.03.003>.

[2]: Bristow, K.L., and G.S. Campbell, 1984. On the relationship between
    incoming solar radiation and daily maximum and minimum temperature.
    Agricultural and Forest Meteorology, 31:159-166.

[3]: Running, S.W., R.R. Nemani, and R.D. Hungerford, 1987. Extrapolation of
    synoptic meteorological data in mountainous terrain and its use for
    simulating forest evaporation and photosynthesis. Canadian Journal of
    Forest Research, 17:472-483.

[4]: Glassy, J.M., and S.W. Running, 1994. Validating diurnal climatology of
    the MT-CLIM model across a climatic gradient in Oregon. Ecological
    Applications, 4(2):248-257.

[5]: Kimball, J.S., S.W. Running, and R. Nemani, 1997. An improved method for
    estimating surface humidity from daily minimum temperature. Agricultural
    and Forest Meteorology, 85:87-98.

[6]: Thornton, P.E., and S.W. Running, 1999. An improved algorithm for
    estimating incident daily solar radiation from measurements of
    temperature, humidity, and precipitation. Agricultural and Forest
    Meteorology, 93:211-228.

[7]: Bohn, T. J., K. M. Whitney, G. Mascaro, and E. R. Vivoni, 2019. A
    deterministic approach for approximating the diurnal cycle of
    precipitation for large-scale hydrological simulations. Journal of
    Hydrometeorology (accepted). doi: 10.1175/JHM-D-18-0203.1.


