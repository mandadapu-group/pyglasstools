R""" Thermodynamic observables.
"""
from pyglasstools.thermo import _thermo
from pyglasstools import _pyglasstools, comm, rank, size, solvers_list
import numpy as np

def initialize_global(names,dim):
    list_obs = {}
    if any("virialstress" in s for s in names):
        list_obs['virialstress'] = _thermo.GlobalVirialStress("virialstress", "2-TENSOR", False, dim)
    elif any("kineticstress" in s for s in names):
        list_obs['kineticstress'] = _thermo.GlobalKineticStress("kineticstress", "2-TENSOR", False, dim)
    elif any("borntensor" in s for s in names):
        list_obs['borntensor'] = _thermo.GlobalBornTensor("borntensor", "4-TENSOR", False, dim)#rcut
    return list_obs

class calculator(object):
    global solvers_list

    def __init__(self, sysdata, potential):
        #Initialize system data and pair potential of the system
        self.sysdata = sysdata;
        self.potential = potential;
        self.manager = _pyglasstools.Manager();
        self.thermocalculator = _thermo.ThermoCalculator(sysdata._getParticleSystem(),potential._getPairPotential())#,self.manager,comm)
        solvers_list.append(self)
   
    def add_observables(self, observables):
        for name in observables:
            self.thermocalculator.addObservable(observables[name])
    def run(self):
        self.thermocalculator.compute()
    
    def update(self,frame_num):
        self.sysdata.update(frame_num); #Let's try and move it up? Have it story current frame number . . .
        self.thermocalculator.setSystemData(self.sysdata._getParticleSystem())
