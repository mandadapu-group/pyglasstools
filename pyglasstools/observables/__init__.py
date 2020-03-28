R""" Pair potentials.
"""
from pyglasstools.observables import _observables
import pyglasstools
import numpy as np

class virialstress(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, dim):
        self.Tv = _observables.VirialStress(dim)#rcut
    
    #Redefine attributes so that it directly access SimBox C++ class 
    #attributes
    def __getattr__(self,attr):
            orig_attr = self.Tv.__getattribute__(attr)
            if callable(orig_attr):
                def hooked(*args, **kwargs):
                    self.pre()
                    result = orig_attr(*args, **kwargs)
                    # prevent Tv from becoming unwrapped
                    if result == self.Tv:
                        return self
                    self.post()
                    return result
                return hooked
            else:
                return orig_attr
    
    def _getObservable(self):
        return self.Tv

class kineticstress(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, dim):
        self.Tk = _observables.KineticStress(dim)#rcut
    
    #Redefine attributes so that it directly access SimBox C++ class 
    #attributes
    def __getattr__(self,attr):
            orig_attr = self.Tk.__getattribute__(attr)
            if callable(orig_attr):
                def hooked(*args, **kwargs):
                    self.pre()
                    result = orig_attr(*args, **kwargs)
                    # prevent Tk from becoming unwrapped
                    if result == self.Tk:
                        return self
                    self.post()
                    return result
                return hooked
            else:
                return orig_attr
    
    def _getObservable(self):
        return self.Tk
