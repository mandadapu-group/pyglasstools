R""" Pair potentials.
"""
from pyglasstools.potential import _potential
import pyglasstools
from pyglasstools import globalpotential
import numpy as np

class lj(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, eps, rcut, name="lennard-jones",mode=None):
        
        self.name = "{}+{}".format(name,mode)
        # create the c++ mirror class
        if (mode == "truncated"):
            self.cpppairpotential = _potential.PairPotentialLJ(rcut,[eps]);
        elif (mode == "force-shifted"):
            self.cpppairpotential = _potential.PairPotentialForceShiftedLJ(rcut,[eps]);
        elif (mode == None):
            raise NameError('Please select a Lennard-Jones potential available modes: truncated and force-shifted are available')
        else:
            raise NameError('Lennard jones potential mode not recognized. Only: truncated and force-shifted are available')
        
        pyglasstools.set_potential(self)

    def get_potentialname(self):
        return "lennard-jones"
    
    def set_diameters(self,diameter_i, diameter_j):
        self.cpppairpotential.di = diameter_i
        self.cpppairpotential.dj = diameter_j
    
    def get_diameters(self):
        return [self.cpppairpotential.di, self.cpppairpotential.dj]
    
    def set_rij(self,r_ij):
        self.cpppairpotential.rij = r_ij.astype('float64')
    
    def get_rij(self):
        return self.cpppairpotential.rij
    
    def set_scaledrcut(self, rcut):
        self.cpppairpotential.scaled_rcut = rcut
    def get_scaledrcut(self):
        return self.cpppairpotential.scaled_rcut

    def set_eps(self, eps):
        self.cpppairpotential.params = [eps]
    def get_eps(self):
        return self.cpppairpotential.params[0]
    
    def get_pairforce(self):
        return self.cpppairpotential.getPairForce()

class polydisperse12(object):
    def __init__(self, v0 =1.0, eps=0.2, rcut=1.25, name="polydisperse-12"):
        self.name = name
        c0 =  -28.0*v0/rcut**12;
        c1 =  48.0*v0/rcut**14;
        c2 =  -21.0*v0/rcut**16;
        self.cpppairpotential = _potential.PairPotentialPoly12(rcut,[v0,eps,c0,c1,c2])
        pyglasstools.set_potential(self)
        
    def _getPairPotential(self):
        return self.cpppairpotential
    def get_potentialname(self):
        return self.name
    def set_diameters(self,diameter_i, diameter_j):
        self.cpppairpotential.di = diameter_i
        self.cpppairpotential.di = diameter_j
    
    def get_diameters(self):
        return [self.cpppairpotential.di, self.cpppairpotential.dj]
    
    def set_rij(self,r_ij):
        self.cpppairpotential.rij = r_ij.astype('float64')
    
    def get_rij(self):
        return self.cpppairpotential.rij
    def set_scaledrcut(self, rcut):
        self.cpppairpotential.scaledrcut = rcut
    def get_scaledrcut(self):
        return self.cpppairpotential.scaledrcut
    
    def set_v0(self, v0):
        self.cpppairpotential.params = [v0,self.cpppairpotential.param[1]]
    def get_v0(self):
        return self.cpppairpotential.params[0]
    
    def set_eps(self, eps):
        self.cpppairpotential.params = [self.cpppairpotential.param[0],eps]
    def get_eps(self):
        return self.cpppairpotential.params[1]
    
    def get_pairforce(self):
        return self.cpppairpotential.getPairForce()

class polydisperse18(object):
    def __init__(self, v0 =1.0, eps=0.0, rcut=1.25, name="polydisperse-18"):
        self.name = name
        c0 =  -55.0*v0/rcut**18;
        c1 =  99.0*v0/rcut**20;
        c2 =  -45.0*v0/rcut**22;
        self.cpppairpotential = _potential.PairPotentialPoly18(rcut,[v0,eps,c0,c1,c2])
        pyglasstools.set_potential(self)
        
    def _getPairPotential(self):
        return self.cpppairpotential
    def get_potentialname(self):
        return self.name
    def set_diameters(self,diameter_i, diameter_j):
        self.cpppairpotential.di = diameter_i
        self.cpppairpotential.di = diameter_j
    
    def get_diameters(self):
        return [self.cpppairpotential.di, self.cpppairpotential.dj]
    
    def set_rij(self,r_ij):
        self.cpppairpotential.rij = r_ij.astype('float64')
    
    def get_rij(self):
        return self.cpppairpotential.rij
    def set_scaledrcut(self, rcut):
        self.cpppairpotential.scaledrcut = rcut
    def get_scaledrcut(self):
        return self.cpppairpotential.scaledrcut
    
    def set_v0(self, v0):
        self.cpppairpotential.params = [v0,self.cpppairpotential.param[1]]
    def get_v0(self):
        return self.cpppairpotential.params[0]
    
    def set_eps(self, eps):
        self.cpppairpotential.params = [self.cpppairpotential.param[0],eps]
    def get_eps(self):
        return self.cpppairpotential.params[1]
    
    def get_pairforce(self):
        return self.cpppairpotential.getPairForce()

class polydisperselj(object):
    def __init__(self, v0 = 1.0, eps=0.2, rcut=2.5, name="polydisperse-lj"):
        self.name = name
        c0 =  -28.0*v0/rcut**12+10.0/rcut**6;
        c1 =  48.0*v0/rcut**14-15.0/rcut**8;
        c2 =  -21.0*v0/rcut**16+6.0/rcut**10;
        self.cpppairpotential = _potential.PairPotentialPolyLJ(rcut,[v0,eps,c0,c1,c2])
        pyglasstools.set_potential(self)
        
    def _getPairPotential(self):
        return self.cpppairpotential
    def get_potentialname(self):
        return self.name
    def set_diameters(self,diameter_i, diameter_j):
        self.cpppairpotential.di = diameter_i
        self.cpppairpotential.di = diameter_j
    
    def get_diameters(self):
        return [self.cpppairpotential.di, self.cpppairpotential.dj]
    
    def set_rij(self,r_ij):
        self.cpppairpotential.rij = r_ij.astype('float64')
    
    def get_rij(self):
        return self.cpppairpotential.rij
    def set_scaledrcut(self, rcut):
        self.cpppairpotential.scaledrcut = rcut
    def get_scaledrcut(self):
        return self.cpppairpotential.scaledrcut
    
    def set_v0(self, v0):
        self.cpppairpotential.params = [v0,self.cpppairpotential.param[1]]
    def get_v0(self):
        return self.cpppairpotential.params[0]
    
    def set_eps(self, eps):
        self.cpppairpotential.params = [self.cpppairpotential.param[0],eps]
    def get_eps(self):
        return self.cpppairpotential.params[1]
    
    def get_pairforce(self):
        return self.cpppairpotential.getPairForce()
