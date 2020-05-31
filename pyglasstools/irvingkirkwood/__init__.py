R""" Irving-Kirkwood coarse-grained fields.
"""
from pyglasstools.irvingkirkwood import _irvingkirkwood
from pyglasstools import _pyglasstools, comm, rank, size, solvers_list
import numpy as np

def initialize_field(names,dim,gridpoints):
    list_obs = {}
    if any("virialstress" in s for s in names):
        if dim == 2:
            list_obs['virialstress'] = _irvingkirkwood.VirialStressField2D("virialstress", "2-TENSOR", False,len(gridpoints))
            list_obs['virialstress'].addGridpoints(gridpoints)
        elif dim == 3:
            list_obs['virialstress'] = _irvingkirkwood.VirialStressField3D("virialstress", "2-TENSOR", False,len(gridpoints))
            list_obs['virialstress'].addGridpoints(gridpoints)
    if any("kineticstress" in s for s in names):
        if dim == 2:
            list_obs['kineticstress'] = _irvingkirkwood.KineticStressField2D("kineticstress", "2-TENSOR", True,len(gridpoints))
            list_obs['kineticstress'].addGridpoints(gridpoints)
        elif dim == 3:
            list_obs['kineticstress'] = _irvingkirkwood.KineticStressField3D("kineticstress", "2-TENSOR", True,len(gridpoints))
            list_obs['kineticstress'].addGridpoints(gridpoints)
    return list_obs

class ikcalculator(object):
    global solvers_list

    def __init__(self, sysdata, potential, cgfunc, dx):
        #Initialize system data and pair potential of the system
        self.sysdata = sysdata;
        self.potential = potential;
        self.cgfunc = cgfunc
        self.manager = _pyglasstools.Manager();
        
        nmax = int(sysdata.simbox.boxsize[0]/dx)
        points = np.linspace(-sysdata.simbox.boxsize[0]/2.0,+sysdata.simbox.boxsize[0]/2.0,nmax)
        self.gridpoints = []
        for i in range(len(points)): 
            for j in range(len(points)):
                self.gridpoints.append(np.array([points[i],points[j],0]).astype(np.float64))
        self.gridpoints = np.asarray(self.gridpoints,dtype=np.float64)
        self.gridsize = len(self.gridpoints) 
        
        self.calculator = _irvingkirkwood.IrvingKirkwood(sysdata._getParticleSystem(),potential._getPairPotential(),cgfunc._getCGFunc(),self.gridpoints)
        solvers_list.append(self)
   
    def add_observables(self, observables):
        for name in observables:
            self.calculator.addObservable(observables[name])
    def run(self):
        self.calculator.compute()
    
    def update(self,frame_num):
        self.sysdata.update(frame_num); #Let's try and move it up? Have it story current frame number . . .
        self.calculator.setSystemData(self.sysdata._getParticleSystem())
