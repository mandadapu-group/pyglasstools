R""" Thermodynamic observables.
"""
from pyglasstools.thermo import _thermo
import pyglasstools
from pyglasstools import _pyglasstools, comm, rank, size, solvers_list
import numpy as np

def initialize_global(names,dim):
    list_obs = {}
    if any("virialstress" in s for s in names):
        if any("elastic" in s for s in names) is True:
            if dim == 2:
                list_obs['elasticvirialstress'] = _thermo.GlobalElasticVirialStress2D("elasticvirialstress", "2-TENSOR", False)
            elif dim == 3:
                list_obs['elasticvirialstress'] = _thermo.GlobalElasticVirialStress3D("elasticvirialstress", "2-TENSOR", False)
        if any("elastic" in s for s in names) is False:
            if dim == 2:
                list_obs['virialstress'] = _thermo.GlobalVirialStress2D("virialstress", "2-TENSOR", False)
            elif dim == 3:
                list_obs['virialstress'] = _thermo.GlobalVirialStress3D("virialstress", "2-TENSOR", False)
    if any("kineticstress" in s for s in names):
        if dim == 2:
            list_obs['kineticstress'] = _thermo.GlobalKineticStress2D("kineticstress", "2-TENSOR", False)
        elif dim == 3:
            list_obs['kineticstress'] = _thermo.GlobalKineticStress3D("kineticstress", "2-TENSOR", False)
    if any("borntensor" in s for s in names):
        if dim == 2:
            list_obs['borntensor'] = _thermo.GlobalBornTensor2D("borntensor", "4-TENSOR", False)
        elif dim == 3:
            list_obs['borntensor'] = _thermo.GlobalBornTensor3D("borntensor", "4-TENSOR", False)
    return list_obs

class calculator(object):
    global solvers_list

    def __init__(self):
        #Initialize system data and pair potential of the system
        self.manager = _pyglasstools.Manager();
        self.thermocalculator = _thermo.ThermoCalculator(   pyglasstools.get_sysdata().cppparticledata,
                                                            pyglasstools.get_potential().cpppairpotential)
        solvers_list.append(self)
   
    def add_observables(self, observables):
        for name in observables:
            self.thermocalculator.addObservable(observables[name])
    def run(self):
        self.thermocalculator.compute()
    
    def update(self,frame_num):
        pyglasstools.update_sysdata(frame_num); #Let's try and move it up? Have it story current frame number . . .
        self.thermocalculator.setSystemData(pyglasstools.get_sysdata().cppparticledata)
