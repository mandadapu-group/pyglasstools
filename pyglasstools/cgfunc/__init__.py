R""" Pair potentials.
"""
from pyglasstools.cgfunc import _cgfunc
import pyglasstools
import numpy as np

class octic(object):
    R""" A truncated 8-th order polynomial coarse grain function.

    """
    def __init__(self, rcut, mode=None):

        self.rcut = rcut
        # create the c++ mirror class
        self.cgfunc = _cgfunc.CGFuncOctic(self.rcut);
    def _getCGFunc(self):
        return self.cgfunc

    def set_rcut(self, rcut):
        self.rcut = rcut
        self.cgfunc.setRcut(rcut)
    def get_rcut(self):
        return self.cgfunc.getRcut()
    
    def set_x(self,x):
        self.cgfunc.setX(x.astype('float64'))
    def get_x(self):
        return self.cgfunc.getX()
    
    def set_ri(self,ri):
        self.cgfunc.setRi(ri.astype('float64'))
    def get_ri(self):
        return self.cgfunc.getRi()
    
    def set_dr(self,dr):
        self.cgfunc.setRij(dr.astype('float64'))
    def get_dr(self):
        return self.cgfunc.getRij()
    
    def get_deltafunc(self):
        return self.cgfunc.getDeltaFunc()
    def get_objfunc(self,s):
        return self.cgfunc.getObjFunc(s)
    def get_bondfunc(self):
        return self.cgfunc.getBondFunc()
