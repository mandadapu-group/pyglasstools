R""" Pair potentials.
"""
from pyglasstools.observables import _observables
import pyglasstools
import numpy as np

def initialize_global(names):
    list_obs = {}
    if any("g_virialstress" in s for s in names):
        list_obs['g_virialstress'] = virialstress(dim=2)
    elif any("g_kineticstress" in s for s in names):
        list_obs['g_kineticstress'] = virialstress(dim=2)
    return list_obs

def initialize_field(names,length):
    list_obs = {}
    if any("f_virialstress" in s for s in names):
        list_obs['f_virialstress'] = virialstressfield(2,length)
    elif any("f_kineticstress" in s for s in names):
        list_obs['f_kineticstress'] = kineticstressfield(2,length)
    return list_obs

class virialstress(object):
    def __init__(self, dim):
        self.Tv = _observables.GlobalVirialStress("g_virialstress", "2-TENSOR", False, True, dim)#rcut
    def _getObservable(self):
        return self.Tv
    def getVal(self):
        return self.Tv.val
class kineticstress(object):
    def __init__(self, dim):
        self.Tk = _observables.GlobalKineticStress("g_kineticstress", "2-TENSOR", True, False, dim)#rcut
    
    def _getObservable(self):
        return self.Tk
    def getVal(self):
        return self.Tk.val

class bornstiffness(object):
    def __init__(self, dim):
        self.CB = _observables.GlobalBornTensor("g_borntensor", "4-TENSOR", False, True, dim)#rcut
    
    #Redefine attributes so that it directly access SimBox C++ class 
    
    def _getObservable(self):
        return self.CB

class virialstressfield(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, dim, gridsize):
        self.Tv = _observables.VirialStressField("f_virialstress", "2-TENSOR", False, True, dim,gridsize)
    def _getObservable(self):
        return self.Tv
    def getVal(self):
        return self.Tv.val

class kineticstressfield(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, dim,gridsize):
        self.Tk = _observables.KineticStressField("f_kineticstress", "2-TENSOR", True, False, dim,gridsize)#rcut
    
    def _getObservable(self):
        return self.Tk
    def getVal(self):
        return self.Tk.val

class densityfield(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, dim, gridsize):
        self.Tv = _observables.VirialStressField("Virial Stress Field", "2-TENSOR", False, True, dim,gridsize)
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
