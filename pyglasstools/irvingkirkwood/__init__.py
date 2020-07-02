R""" Irving-Kirkwood coarse-grained fields.
"""
from pyglasstools.irvingkirkwood import _irvingkirkwood
from pyglasstools import _pyglasstools, comm, rank, size, solvers_list
import numpy as np

def initialize_field(names,dim,gridpoints):
    list_obs = {}
    if any("virialstress" in s for s in names):
        if any("elastic" in s for s in names):
            if dim == 2:
                list_obs['elasticvirialstress'] = _irvingkirkwood.ElasticVirialStressField2D("elasticvirialstress", "2-TENSOR", False,len(gridpoints))
                list_obs['elasticvirialstress'].addGridpoints(gridpoints)
            elif dim == 3:
                list_obs['elasticvirialstress'] = _irvingkirkwood.ElasticVirialStressField3D("elasticvirialstress", "2-TENSOR", False,len(gridpoints))
                list_obs['elasticvirialstress'].addGridpoints(gridpoints)
        if all("elastic" in s for s in names) is False:
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
        
        #Maybe move it to C++ side . . . .
        globalsize = 0
        if rank == 0:
            nmax = int(sysdata.simbox.boxsize[0]/dx)
            points = np.linspace(-sysdata.simbox.boxsize[0]/2.0,+sysdata.simbox.boxsize[0]/2.0,nmax)
            gridpoints = []
            for i in range(len(points)): 
                for j in range(len(points)):
                    gridpoints.append(np.asarray([points[i],points[j],0],dtype=np.float64))
            globalsize = len(gridpoints)
            self.gridpoints = np.array_split(np.asarray(gridpoints,dtype=np.float64),size)
            del points, nmax
        else:
            self.gridpoints = np.array_split([np.zeros(size)],size) 

        #Scatter and reshape
        self.gridpoints = comm.scatter_v(self.gridpoints,0)
        self.gridpoints = np.reshape(self.gridpoints,(len(self.gridpoints),3))
        
        #Broadcast the true size
        self.gridsize =  comm.bcast(globalsize,0); 
        self.calculator = _irvingkirkwood.IrvingKirkwood(sysdata._getParticleSystem(),potential._getPairPotential(),cgfunc._getCGFunc(),comm)
        solvers_list.append(self)
   
    def add_observables(self, observables):
        for name in observables:
            self.calculator.addObservable(observables[name])
    def run(self):
        self.calculator.compute(self.gridpoints)
    
    def update(self,frame_num):
        self.sysdata.update(frame_num); #Let's try and move it up? Have it story current frame number . . .
        self.calculator.setSystemData(self.sysdata._getParticleSystem())


class radialcalculator(object):
    global solvers_list

    def __init__(self, sysdata, potential, cgfunc, center=np.array([0,0,0]),spacingtype="normal",rmax=10,rmin=0,dr=0.1,dlnr=0.05):
        #Initialize system data and pair potential of the system
        self.sysdata = sysdata;
        self.potential = potential;
        self.cgfunc = cgfunc
        self.manager = _pyglasstools.Manager();

        #TO DO: Move all of this hassle to C++ side . . . .
        globalsize = 0
        self.center = center
        if rank == 0:
            gridpoints = []
            if spacingtype == "normal":
                r = np.linspace(0,rmax,int(rmax/dr))[1:]
            elif spacingtype == "log":
                r = np.exp(np.linspace(np.log(rmin),np.log(rmax),int((np.log(rmax)-np.log(rmin))/dlnr)))
            
            gridid = [] #id based on radius
            globalr = []  
            globaltheta = []
            gridpoints = []            
            
            #Append the origin
            gridid.append(0)
            globalr.append(0)
            globaltheta.append(0)
            gridpoints.append(self.center)

            for i in range(len(r)): 
                dtheta = dr/r[i]
                theta = np.linspace(0,2*np.pi,int(2*np.pi/dtheta))
                x = r[i]*np.cos(theta) 
                y = r[i]*np.sin(theta)
                for j in range(len(theta)):
                    gridid.append(i+1)
                    globalr.append(r[i])
                    globaltheta.append(theta[j])
                    gridpoints.append(np.asarray([x[j],y[j],0],dtype=np.float64)+self.center)
            
            globalsize = len(gridpoints)
            #print(np.shape(gridpoints),np.shape(center))
            self.gridid = np.array_split(np.asarray(gridid,dtype=np.float64),size)
            self.globalr = np.array_split(np.asarray(globalr,dtype=np.float64),size)
            self.globaltheta = np.array_split(np.asarray(globaltheta,dtype=np.float64),size)
            self.gridpoints = np.array_split(np.asarray(gridpoints,dtype=np.float64),size)
            del r, theta, dtheta, x, y, gridpoints, gridid, globalr, globaltheta
        else:
            self.gridid = np.array_split(np.zeros(size),size) 
            self.globalr = np.array_split(np.zeros(size),size) 
            self.globaltheta = np.array_split(np.zeros(size),size) 
            self.gridpoints = np.array_split([np.zeros(size)],size) 

        #Scatter and reshape
        self.gridid = comm.scatter_v(self.gridid,0)
        self.gridid = np.reshape(self.gridid,(len(self.gridid),))
        #Scatter and reshape
        self.globalr = comm.scatter_v(self.globalr,0)
        self.globalr = np.reshape(self.globalr,(len(self.globalr),))
        #Scatter and reshape
        self.globaltheta = comm.scatter_v(self.globaltheta,0)
        self.globaltheta = np.reshape(self.globaltheta,(len(self.globaltheta),))
        #Scatter and reshape
        self.gridpoints = comm.scatter_v(self.gridpoints,0)
        self.gridpoints = np.reshape(self.gridpoints,(len(self.gridpoints),3))
        
        #Broadcast the true size
        self.gridsize =  comm.bcast(globalsize,0); #del globalsize
        self.calculator = _irvingkirkwood.IrvingKirkwood(sysdata._getParticleSystem(),potential._getPairPotential(),cgfunc._getCGFunc(),comm)
        solvers_list.append(self)

    def set_center(self,center):
        self.center = center
        self.gridpoints[0:] += center
    
    def add_observables(self, observables):
        for name in observables:
            self.calculator.addObservable(observables[name])
    
    def run(self):
        self.calculator.compute(self.gridpoints)
    
    def update(self,frame_num):
        self.sysdata.update(frame_num); #Let's try and move it up? Have it story current frame number . . .
        self.calculator.setSystemData(self.sysdata._getParticleSystem())
