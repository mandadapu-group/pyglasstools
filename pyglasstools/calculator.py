from pyglasstools import _pyglasstools
import numpy as np

class ikglobal(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, sysdata,potential, vol):
        self.irvingkirkwood = _pyglasstools.GlobalCalculator(sysdata._getParticleSystem(),potential._getPairPotential());
        self.vol = vol
    def add_observable(self,obs):
        self.irvingkirkwood.addObservable(obs._getObservable());
    def compute(self):
        self.irvingkirkwood.compute();
    def get_virialstressvalue(self):
        return self.irvingkirkwood.getGlobalObservable("Virial Stress")/self.vol;
    def get_kineticstressvalue(self):
        return self.irvingkirkwood.getGlobalObservable("Kinetic Stress")/self.vol;

class iklocal(object):
    def __init__(self, sysdata, potential, cgfunc):
        self.irvingkirkwood = _pyglasstools.LocalCalculator(sysdata._getParticleSystem(),potential._getPairPotential(), cgfunc._getCGFunc());
    def add_observable(self,obs):
        self.irvingkirkwood.addObservable(obs._getObservable());
    def compute(self,gridpoints):
        self.irvingkirkwood.computelocal(gridpoints);
    def get_virialstress(self):
        return np.array(self.irvingkirkwood.getField("Virial Stress Field")).astype('float64')
    def get_kineticstress(self):
        return np.array(self.irvingkirkwood.getField("Kinetic Stress Field")).astype('float64')
    def get_density(self):
        return np.array(self.irvingkirkwood.getField("Density Field")).astype('float64')
