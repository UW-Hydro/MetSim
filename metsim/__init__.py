"""
MetSim global module
"""

from .metsim import MetSim

from . import _version
__version__ = _version.get_versions()['version']
