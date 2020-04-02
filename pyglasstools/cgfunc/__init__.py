R""" Pair potentials.
"""
from pyglasstools.cgfunc import _cgfunc
import pyglasstools
import numpy as np

class octic(object):
    R""" A truncated 8-th order polynomial coarse grain function.

    """
    def __init__(self, rcut, mode=None):

        # create the c++ mirror class
        self.cgfunc = _cgfunc.CGFuncOctic(rcut);
    
    def _getCGFunc(self):
        return self.cgfunc

    def set_rcut(self, rcut):
        self.cgfunc.cg_rcut = rcut
    def get_rcut(self):
        return self.cgfunc.cg_rcut

    def set_x(self,x):
        self.cgfunc.x = x.astype('float64')
    def get_x(self):
        return self.cgfunc.x
    
    def set_ri(self,ri):
        self.cgfunc.ri = ri
    def get_ri(self):
        return self.cgfunc.ri

    def set_rij(self,rij):
        self.cgfunc.rij = rij.astype('float64')
    def get_rij(self):
        return self.cgfunc.rij
    
    def get_deltafunc(self):
        return self.cgfunc.getDeltaFunc()
    def get_objfunc(self,s):
        return self.cgfunc.getObjFunc(s)
    def get_bondfunc(self):
        return self.cgfunc.getBondFunc()
