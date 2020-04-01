R""" Pair potentials.
"""
from pyglasstools.observables import _observables
import pyglasstools
import numpy as np

class virialstress(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, dim):
        self.Tv = _observables.GlobalVirialStress("Virial Stress", "TENSOR", False, True, dim)#rcut
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
        self.Tk = _observables.GlobalKineticStress("Kinetic Stress", "TENSOR", True, False, dim)#rcut
    
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

class virialstressfield(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, dim, gridsize):
        self.Tv = _observables.VirialStressField("Virial Stress Field", "TENSOR", False, True, dim,gridsize)
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

class kineticstressfield(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, dim,gridsize):
        self.Tk = _observables.KineticStressField("Kinetic Stress Field", "TENSOR", True, False, dim,gridsize)#rcut
    
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

class densityfield(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, dim, gridsize):
        self.Tv = _observables.VirialStressField("Virial Stress Field", "TENSOR", False, True, dim,gridsize)
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

class kineticstressfield(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, dim,gridsize):
        self.rho = _observables.DensityField("Density Field", "SCALAR", True, False, dim, gridsize)
    
    #Redefine attributes so that it directly access SimBox C++ class 
    #attributes
    def __getattr__(self,attr):
            orig_attr = self.rho.__getattribute__(attr)
            if callable(orig_attr):
                def hooked(*args, **kwargs):
                    self.pre()
                    result = orig_attr(*args, **kwargs)
                    # prevent rho from becoming unwrapped
                    if result == self.rho:
                        return self
                    self.post()
                    return result
                return hooked
            else:
                return orig_attr
    
    def _getObservable(self):
        return self.rho
