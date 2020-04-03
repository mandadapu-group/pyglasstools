R""" Pair potentials.
"""
from pyglasstools.potential import _potential
import numpy as np

class lj(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, eps, rcut, name="lennard-jones",mode=None):
        
        self.name = "{}+{}".format(name,mode)
        # create the c++ mirror class
        if (mode == "truncated"):
            self.lennardjones = _potential.PairPotentialLJ(rcut,[eps]);
        elif (mode == "force-shifted"):
            self.lennardjones = _potential.PairPotentialForceShiftedLJ(rcut,[eps]);
        elif (mode == None):
            raise NameError('Please select a Lennard-Jones potential available modes: truncated and force-shifted are available')
        else:
            raise NameError('Lennard jones potential mode not recognized. Only: truncated and force-shifted are available')
    
    def _getPairPotential(self):
        return self.lennardjones
    def get_potentialname(self):
        return "lennard-jones"
    
    def set_diameters(self,diameter_i, diameter_j):
        self.lennardjones.di = diameter_i
        self.lennardjones.dj = diameter_j
    
    def get_diameters(self):
        return [self.lennardjones.di, self.lennardjones.dj]
    
    def set_rij(self,r_ij):
        self.lennardjones.rij = r_ij.astype('float64')
    
    def get_rij(self):
        return self.lennardjones.rij
    
    def set_scaledrcut(self, rcut):
        self.lennardjones.scaled_rcut = rcut
    def get_scaledrcut(self):
        return self.lennardjones.scaled_rcut

    def set_eps(self, eps):
        self.lennardjones.params = [eps]
    def get_eps(self):
        return self.lennardjones.params[0]
    
    def get_pairforce(self):
        return self.lennardjones.getPairForce()

class polydisperse12(object):
    def __init__(self, v0 =1.0, eps=0.2, rcut=1.25, name="polydisperse-12"):
        self.name = name
        c0 =  -28.0*v0/rcut**12;
        c1 =  48.0*v0/rcut**14;
        c2 =  -21.0*v0/rcut**16;
        self.polydisperse = _potential.PairPotentialPoly12(rcut,[v0,eps,c0,c1,c2])
    
    def _getPairPotential(self):
        return self.polydisperse
    def get_potentialname(self):
        return self.name
    def set_diameters(self,diameter_i, diameter_j):
        self.polydisperse.di = diameter_i
        self.polydisperse.di = diameter_j
    
    def get_diameters(self):
        return [self.polydisperse.di, self.polydisperse.dj]
    
    def set_rij(self,r_ij):
        self.polydisperse.rij = r_ij.astype('float64')
    
    def get_rij(self):
        return self.polydisperse.rij
    def set_scaledrcut(self, rcut):
        self.polydisperse.scaledrcut = rcut
    def get_scaledrcut(self):
        return self.polydisperse.scaledrcut
    
    def set_v0(self, v0):
        self.polydisperse.params = [v0,self.polydisperse.param[1]]
    def get_v0(self):
        return self.polydisperse.params[0]
    
    def set_eps(self, eps):
        self.polydisperse.params = [self.polydisperse.param[0],eps]
    def get_eps(self):
        return self.polydisperse.params[1]
    
    def get_pairforce(self):
        return self.polydisperse.getPairForce()
