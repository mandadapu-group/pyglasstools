R""" Pair potentials.
"""
from pyglasstools import _pyglasstools
import pyglasstools
import numpy as np

class lj(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, eps, rcut, mode=None):

        self.rcut = rcut
        self.params = [];
        self.params.append(eps)
        self.params = np.array(self.params).astype('float64')
        
        # create the c++ mirror class
        if (mode == "truncated"):
            self.pairpotential = _pyglasstools.PairPotentialLJ(self.rcut,self.params);
        elif (mode == "force-shifted"):
            self.pairpotential = _pyglasstools.PairPotentialForceShiftedLJ(self.rcut,self.params);
        elif (mode == None):
            raise NameError('Please select a Lennard-Jones potential available modes: truncated and force-shifted are available')
        else:
            raise NameError('Lennard jones potential mode not recognized. Only: truncated and force-shifted are available')
    def _getPairPotential(self):
        return self.pairpotential
    def get_potentialname(seld):
        return "lennard-jones"
    def set_diameters(self,diameter_i, diameter_j):
        self.pairpotential.setDiameters(diameter_i,diameter_j)
    
    def get_diameters(self):
        return self.pairpotential.getDiameters()
    
    def set_rij(self,r_ij):
        self.pairpotential.setRij(r_ij.astype('float64'))
    
    def get_rij(self):
        return self.pairpotential.getRij()
    
    def set_rcut(self, rcut):
        self.rcut = rcut
        self.pairpotential.setRcut(rcut)
    def get_rcut(self):
        return self.pairpotential.getRcut()
    
    def set_eps(self, eps):
        self.params[0] = eps
        self.pairpotential.setParams(self.params)
    def get_eps(self):
        return self.pairpotential.getParams()
    
    def get_pairforce(self):
        return self.pairpotential.getPairForce()
