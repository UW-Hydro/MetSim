"""
MetSim global module
"""

from metsim.methods import mtclim

METHODS = {
           'mtclim': mtclim
          }


def multi(dispatch):
    def _inner(arg):
        return _inner.__multi__.get(dispatch(arg), _inner.__multi_error__)(arg)

    _inner.__multi__ = {}
    _inner.__multi_error__ = lambda x : None 
    return _inner


def method(dispatch, key):
    def decorate(fn):
        if key is None:
            dispatch.__multi_error__ = fn
        else:
            dispatch.__multi__[key] = fn
        return dispatch
    return decorate


