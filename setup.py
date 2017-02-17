#!/usr/bin/env python

try:
    from setuptools import setup
except:
    from distutils.core import setup

setup(name='metsim',
      version='0.0.0',
      description='Meteorology Simulator',
      author='Andrew Bennett',
      author_email='bennett.andr@gmail.com',
      packages=['metsim'],
      scripts=['scripts/ms'],
      install_requires=['xarray', 'numba'],
      tests_require=['pytest'],)

