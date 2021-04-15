R""" Pair potentials.
"""
from pyglasstools.irvingkirkwood.cgfunc import _cgfunc
import pyglasstools
import numpy as np


class basecgfunc(object):
    R""" A truncated 8-th order polynomial coarse grain function.

    """
    def _getCGFunc(self):
        return self.cgfunc

    def set_rcut(self, rcut):
        self.cgfunc.setRcut(rcut)
    def get_rcut(self):
        return self.cgfunc.getRcut()

    def get_deltafunc(self,x,ri):
        return self.cgfunc.getDeltaFunc(x,ri)
    def get_objfunc(self,s,x,ri,rij):
        return self.cgfunc.getObjFunc(s,x,ri,rij)
    def get_bondfunc(self,x,ri,rij):
        return self.cgfunc.getBondFunc(x,ri,rij)

class octic(basecgfunc):
    R""" A truncated 8-th order polynomial coarse grain function.

    """
    def __init__(self, rcut=2.5, mode=None,**kwargs):
        # create the c++ mirror class
        if (mode == "fixedpoint"):
            if any("order" in s for s in kwargs):
                order = kwargs["order"]
            self.cgfunc = _cgfunc.CGFuncOcticFixedPoint(order,rcut);
        elif (mode == "adaptive"): 
            if any("relerr" in s for s in kwargs):
                relerr = kwargs["relerr"]
            if any("abserr" in s for s in kwargs):
                abserr = kwargs["abserr"]
            if any("maxiter" in s for s in kwargs):
                maxiter = kwargs["maxiter"]
            self.cgfunc = _cgfunc.CGFuncOcticAdaptive(maxiter,relerr,abserr,rcut);

class rect(basecgfunc):
    R""" A truncated radial heaviside coarse grain function.

    """
    def __init__(self, rcut=2.5, mode=None,**kwargs):
        # create the c++ mirror class
        if (mode == "fixedpoint"):
            if any("order" in s for s in kwargs):
                order = kwargs["order"]
            self.cgfunc = _cgfunc.CGFuncRectFixedPoint(order,rcut);
        elif (mode == "adaptive"): 
            if any("relerr" in s for s in kwargs):
                relerr = kwargs["relerr"]
            if any("abserr" in s for s in kwargs):
                abserr = kwargs["abserr"]
            if any("maxiter" in s for s in kwargs):
                maxiter = kwargs["maxiter"]
            self.cgfunc = _cgfunc.CGFuncRectAdaptive(maxiter,relerr,abserr,rcut);

class mollifier(basecgfunc):
    R""" THe bump/mollifier function

    """
    def __init__(self, rcut=2.5, mode=None,**kwargs):
        # create the c++ mirror class
        if (mode == "fixedpoint"):
            if any("order" in s for s in kwargs):
                order = kwargs["order"]
            self.cgfunc = _cgfunc.CGFuncMollifierFixedPoint(order,rcut);
        elif (mode == "adaptive"): 
            if any("relerr" in s for s in kwargs):
                relerr = kwargs["relerr"]
            if any("abserr" in s for s in kwargs):
                abserr = kwargs["abserr"]
            if any("maxiter" in s for s in kwargs):
                maxiter = kwargs["maxiter"]
            self.cgfunc = _cgfunc.CGFuncMollifierAdaptive(maxiter,relerr,abserr,rcut);
