.. _installing:

Installation
============

Required dependencies
---------------------

- Python 3.6
- `xarray <http://xarray.pydata.org/>`__ (0.9.1 or later)
- `pandas <http://pandas.pydata.org/>`__ (0.19.0 or later)
- `numba <http://numba.pydata.org/>`__ (0.31.0 or later)


Optional dependencies
---------------------

For netCDF and IO
~~~~~~~~~~~~~~~~~

- `netCDF4 <https://github.com/Unidata/netcdf4-python>`__: recommended if you
  want to use xarray for reading or writing netCDF files
- `scipy <http://scipy.org/>`__: used as a fallback for reading/writing netCDF3
- `pydap <http://www.pydap.org/>`__: used as a fallback for accessing OPeNDAP
- `h5netcdf <https://github.com/shoyer/h5netcdf>`__: an alternative library for
  reading and writing netCDF4 files that does not use the netCDF-C libraries
- `pynio <https://www.pyngl.ucar.edu/Nio.shtml>`__: for reading GRIB and other
  geoscience specific file formats

For accelerating xarray and metsim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `bottleneck <https://github.com/kwgoodman/bottleneck>`__: speeds up
  NaN-skipping and rolling window aggregations by a large factor
  (1.0 or later)
- `cyordereddict <https://github.com/shoyer/cyordereddict>`__: speeds up most
  internal operations with xarray data structures (for python versions < 3.5)

For parallel computing
~~~~~~~~~~~~~~~~~~~~~~

- `dask.array <http://dask.pydata.org>`__ (0.9.0 or later): required for
  :ref:`dask`.

For plotting
~~~~~~~~~~~~

Instructions
------------

`metsim` itself is a pure Python package, but its dependencies are not. The
easiest way to get everything installed is to use conda_. To install metsim
with its recommended dependencies using the conda command line tool::

    $ conda install metsim

.. _conda: http://conda.io/

We recommend using the community maintained `conda-forge <https://conda-forge.github.io/>`__ channel if you need difficult\-to\-build dependencies such as cartopy or pynio::

    $ conda install -c conda-forge metsim

New releases may also appear in conda-forge before being updated in the default
channel.

If you don't use conda, be sure you have the required dependencies (numba and xarray) installed first. Then, install metsim with pip::

    $ pip install metsim

To run the test suite after installing metsim, install
`py.test <https://pytest.org>`__ (``pip install pytest``) and run
``py.test --verbose``.
