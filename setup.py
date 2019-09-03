#!/usr/bin/env python

import versioneer
from setuptools import setup

setup(name='metsim',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Meteorology Simulator',
      url='https://github.com/UW-Hydro/MetSim',
      author='Andrew Bennett',
      author_email='bennett.andr@gmail.com',
      packages=['metsim', 'metsim.methods', 'metsim.cli'],
      entry_points={
          'console_scripts': ['ms = metsim.cli.ms:main']},
      install_requires=['xarray>=0.11.0', 'numba', 'numpy', 'pandas',
                        'dask', 'distributed', 'toolz', 'netCDF4', 'scipy'],
      keywords=['meteorology', 'disaggregation', 'hydrology',
                'climate', 'mtclim'],
      tests_require=['pytest'],)
