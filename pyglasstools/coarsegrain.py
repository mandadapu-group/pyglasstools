from pyglasstools import _pyglasstools
import pyglasstools

class ikglobal(object):
    R""" Lennard-Jones pair potential.

    """
    def __init__(self, sysdata,vol):
        self.irvingkirkwood = _pyglasstools.GlobalCalculator(sysdata._getParticleSystem());
        self.vol = vol
    def add_observable(self,obs):
        self.irvingkirkwood.addObservable(obs._getObservable());
    def compute(self):
        self.irvingkirkwood.compute();
    def get_virialstressvalue(self):
        return self.irvingkirkwood.getObservable("Virial Stress")/self.vol;
    def get_kineticstressvalue(self):
        return self.irvingkirkwood.getObservable("Kinetic Stress")/self.vol;
