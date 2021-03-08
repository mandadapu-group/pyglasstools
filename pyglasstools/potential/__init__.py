R""" Pair potentials.
"""
from pyglasstools.potential import _potential
import pyglasstools
from pyglasstools import globalpotential
import numpy as np


class pairpotential(object):
    R""" an interface to layout common methods for a generic pair potential

    """
    def __init__(self):
        self.cpppairpotential = None
        self.name = None
        self.dict_params = None

    def get_potentialname(self):
        return self.name
    
    def set_scaledrcut(self, rcut):
        self.cpppairpotential.scaled_rcut = rcut
    def get_scaledrcut(self):
        return self.cpppairpotential.scaled_rcut

    def set_params(self, val, num):
        if isinstance(num,str):
            self.cpppairpotential.setParams(val,self.dict_params[num])
        elif isinstance(num,int):
            self.cpppairpotential.setParams(val,num)
    def get_params(self,num = None):
        if num == None:
            return self.cpppairpotential.getParams()
        elif isinstance(num,str):
            return self.cpppairpotential.getParams()[self.dict_params[num]]
        elif isinstance(num,int):
            return self.cpppairpotential.getParams()[num]
    def print_params(self):
        params = self.cpppairpotential.getParams()
        for key in self.dict_params:
            print("{}: {} ".format(key, params[self.dict_params[key]]))
    def get_pairforce(self,rij,di,dj):
        if len(rij) == 2:
            rij = np.append(np.array(rij),0)
            sigmasq = self.cpppairpotential.getRcut(rij,di,dj)/self.cpppairpotential.scaled_rcut**2
            return self.cpppairpotential.getPairForceDivR(rij,di,dj)*np.linalg.norm(rij)#/sigmasq**(1/2)
        else:
            rij = np.array(rij)
            sigmasq = self.cpppairpotential.getRcut(rij,di,dj)/self.cpppairpotential.scaled_rcut**2
            return self.cpppairpotential.getPairForceDivR(rij,di,dj)*np.linalg.norm(rij)#/sigmasq**(1/2)
    def get_pairenergy(self,rij,di,dj):
        if len(rij) == 2:
            return self.cpppairpotential.getPairEnergy(np.append(np.array(rij),0),di,dj)
        else:
            return self.cpppairpotential.getPairEnergy(np.array(rij),di,dj)


class lj(pairpotential):
    R""" A Lennard-Jones pair potential
         It includes both truncated and force-shifted versions
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
        self.dict_params = {'eps': 0} 
        pyglasstools.set_potential(self)
    

class polydisperse12(pairpotential):
    R""" Poly12 Pair Potential
    """
    def __init__(self, v0 =1.0, eps=0.2, rcut=1.25):
        self.name = "polydisperse-12"
        c0 =  -28.0*v0/rcut**12;
        c1 =  48.0*v0/rcut**14;
        c2 =  -21.0*v0/rcut**16;
        self.cpppairpotential = _potential.PairPotentialPoly12(rcut,[v0,eps,c0,c1,c2])
        self.dict_params = {'v0': 0, 'eps': 1, 'c0': 2, 'c1': 3,'c2': 4}
        pyglasstools.set_potential(self)

class polydisperse18(pairpotential):
    R""" Poly18 Pair Potential
    """
    def __init__(self, v0 =1.0, eps=0.0, rcut=1.25):
        self.name = "polydisperse-18"
        c0 =  -55.0*v0/rcut**18;
        c1 =  99.0*v0/rcut**20;
        c2 =  -45.0*v0/rcut**22;
        self.cpppairpotential = _potential.PairPotentialPoly18(rcut,[v0,eps,c0,c1,c2])
        self.dict_params = {'v0': 0, 'eps': 1, 'c0': 2, 'c1': 3,'c2': 4}
        pyglasstools.set_potential(self)

class polydisperselj(pairpotential):
    R""" PolyLJ Pair Potential
    """
    def __init__(self, v0 = 1.0, eps=0.2, rcut=2.5):
        self.name = "polydisperse-lj"
        c0 =  -28.0*v0/rcut**12+10.0/rcut**6;
        c1 =  48.0*v0/rcut**14-15.0/rcut**8;
        c2 =  -21.0*v0/rcut**16+6.0/rcut**10;
        self.cpppairpotential = _potential.PairPotentialPolyLJ(rcut,[v0,eps,c0,c1,c2])
        self.dict_params = {'v0': 0, 'eps': 1, 'c0': 2, 'c1': 3,'c2': 4}
        pyglasstools.set_potential(self)

class polydisperse10(pairpotential):
    def __init__(self, v0 =1.0, eps=0.0416667, rcut=1.48):
        self.name = "polydisperse-10"
        c0 =  -(56.0)*v0/(rcut**10);
        c1 =  (140.0)*v0/(rcut**12);
        c2 =  -(120.0)*v0/(rcut**14);
        c3 =  (35.0)*v0/(rcut**16);
        self.cpppairpotential = _potential.PairPotentialPoly10(rcut,[v0,eps,c0,c1,c2,c3])
        self.dict_params = {'v0': 0, 'eps': 1, 'c0': 2, 'c1': 3,'c2': 4,'c3': 5}
        pyglasstools.set_potential(self)

class polydisperse106(pairpotential):
    def __init__(self, v0 =1.0, eps=0.1, rcut=3.0):
        self.name = "polydisperse-106"
        c0 = (-21 + 10*rcut**4)*v0/rcut**10
        c1 =  -((5*(-7 + 3*rcut**4)*v0)/rcut**12)
        c2 = (3*(-5 + 2*rcut**4)*v0)/rcut**14
        self.cpppairpotential = _potential.PairPotentialPoly106(rcut,[v0,eps,c0,c1,c2])
        self.dict_params = {'v0': 0, 'eps': 1, 'c0': 2, 'c1': 3,'c2': 4}
        pyglasstools.set_potential(self)

class polydisperseyukawa(pairpotential):
    def __init__(self, v0=10.0, eps=0.0, rcut=4.0, kappa=1.0):
        self.name = "polydisperse-yukawa"
        c0 = -((np.exp(-kappa*rcut)*(15 + 7*kappa*rcut + kappa**2*rcut**2)*v0)/(8*rcut))
        c1 = (np.exp(-kappa*rcut)*(5 + 5*kappa*rcut + kappa**2*rcut**2)*v0)/(4*rcut**3)
        c2 = -((np.exp(-kappa*rcut)*(3 + 3*kappa*rcut + kappa**2*rcut**2)*v0)/(8*rcut**5))
        self.cpppairpotential = _potential.PairPotentialPolyYukawa(rcut,[v0,kappa,eps,c0,c1,c2])
        self.dict_params = {'v0': 0, 'kappa': 1, 'eps': 2, 'c0': 3, 'c1': 4,'c2': 5}
        pyglasstools.set_potential(self)
