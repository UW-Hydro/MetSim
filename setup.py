#!/usr/bin/env python

try:
    from setuptools import setup
except:
    from distutils.core import setup

setup(name='metsim',
      version='1.0.0',
      description='Meteorology Simulator',
      url='https://github.com/UW-Hydro/MetSim',
      download_url='https://github.com/UW-Hydro/MetSim/archive/v1.0.tar.gz',
      author='Andrew Bennett',
      author_email='bennett.andr@gmail.com',
      packages=['metsim', 'metsim.methods'],
      scripts=['scripts/ms', 'scripts/wrap_metsim_over_time.py'],
      install_requires=['xarray', 'numba'],
      keywords=['meteorology', 'disaggregation', 'hydrology',
                'climate', 'mtclim'],
      tests_require=['pytest'],)
