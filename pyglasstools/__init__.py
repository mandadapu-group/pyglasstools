""" PyGlassTools python API
"""
from pyglasstools import _pyglasstools;
from pyglasstools import utils;
from pyglasstools import observables
from os import path
import numpy as np

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
nprocs = size
#The module should save more than one logger, which analyzes various observables
loggers = [];

#https://stackoverflow.com/questions/45680050/cannot-write-to-shared-mpi-file-with-mpi4py
class MPILogFile(object):
    def __init__(self, comm, filename, mode):
        self.file_handle = MPI.File.Open(comm, filename, mode)
        self.file_handle.Set_atomicity(True)
        self.buffer = bytearray

    def write(self, msg):
        b = bytearray()
        b.extend(map(ord, msg))
        self.file_handle.Write_shared(b)

    def close(self):
        self.file_handle.Sync()
        self.file_handle.Close()

#will automatically create the required calculators . . .
class logger(object):
    global loggers 

    def __init__(self, filename=None, names = None, sysdata = None, pair = None):
        
        #Save filename
        self.filename = filename
        
        #First, parse the list of names based on what type of obsercables they are
        self.__glob_obs_names = [s for s in names if "g_" in s ];
        self.global_obs = [] 
        #self.__nonaffine_obs_names = [s for s in names if "na_" in s ];
        #self.nonaffine_obs = [] 
        
        #Initialize system data and pair potential of the system
        self.sysdata = sysdata;
        self.pair = pair;

        #Initialize the "calculators"
        self.ikglobal = None;
        
        if not (not self.__glob_obs_names):
            #Construct the global calculator
            self.ikglobal = _pyglasstools.GlobalCalculator(self.sysdata._getParticleSystem(),self.pair._getPairPotential())
            #Then, add observables
            self.global_obs = observables.initialize_global(self.__glob_obs_names)
            for name in self.global_obs:
                self.ikglobal.addObservable(self.global_obs[name]._getObservable())
        
        #Once initialization is done, the logger adds itself to the global list of available loggers 
        loggers.append(self)

        #Initialize file to save:
        self.file = open(self.filename,"w+")

        #Create column headers
        self.file.write("Frame ")
        for name in self.__glob_obs_names:
            self.file.write("{} ".format(name))
        self.file.write("\n")

    def run(self): 
        if not (not self.__glob_obs_names):
            self.ikglobal.compute()
    
    def save(self,frame_num):
        self.file.write("{:d} ".format(frame_num))
        for name in self.__glob_obs_names:
            val = 0
            if "g_virialstress" in name:
                val = self.global_obs['g_virialstress'].getVal()[int(name[-2]),int(name[-1])]
            elif "g_kineticstress" in name:
                val = self.global_obs['g_kineticstress'].getVal()[int(name[-2]),int(name[-1])]
            elif "g_borntensor" in name:
                val = self.global_obs['g_borntensor'].getVal()[int(name[-2]),int(name[-1])]
            self.file.write("{:.12f} ".format(val/self.sysdata.simbox.vol))
        self.file.write("\n")
    
    def update(self,frame_num):
        self.sysdata.update(frame_num);
        #Update the global calculator
        self.ikglobal = _pyglasstools.GlobalCalculator(self.sysdata._getParticleSystem(),self.pair._getPairPotential())
        #Then, add observables
        self.global_obs = observables.initialize_global(self.__glob_obs_names)
        for name in self.global_obs:
            self.ikglobal.addObservable(self.global_obs[name]._getObservable())

#will automatically create the required calculators . . .
class field_logger(object):
    global loggers 

    def __init__(self, fileprefix=None, names = None, sysdata = None, pair = None, cgfunc = None, dx=1.0):
        
        #Field logger a
        if rank == 0: 
            nmax = int(sysdata.simbox.boxsize[0]/dx)
            points = np.linspace(-sysdata.simbox.boxsize[0]/2.0,+sysdata.simbox.boxsize[0]/2.0,nmax)
            x = []
            for i in range(len(points)): 
                for j in range(len(points)):
                    x.append(np.array([points[i],points[j],0]).astype(np.float64))
            gridpoints = np.asarray(x,dtype=np.float64)
            self.gridsize = len(gridpoints) 
            
            # determine the size of each sub-task
            ave, res = divmod(len(gridpoints), nprocs)
            counts = [ave + 1 if p < res else ave for p in range(nprocs)]


            # determine the starting and ending indices of each sub-task
            starts = [sum(counts[:p]) for p in range(nprocs)]
            ends = [sum(counts[:p+1]) for p in range(nprocs)]

            # converts gridpoints into a list of arrays 
            gridpoints = [gridpoints[starts[p]:ends[p]] for p in range(nprocs)]
        else:
            gridpoints = None
            self.gridsize = None
        self.gridpoints = comm.scatter(gridpoints, root=0)
        self.gridsize = comm.bcast(self.gridsize,root=0) 
        #Delete unnecessarry objects
        del gridpoints
        if rank == 0:
            del nmax,points,ave,res,counts,starts,ends
        
        #Save filename
        self.fileprefix = fileprefix
        
        #First, parse the list of names based on what type of obsercables they are
        self.__field_obs_names = [s for s in names if "f_" in s ];
        self.field_obs = [] 
        
        #Initialize system data and pair potential of the system
        self.sysdata = sysdata;
        self.pair = pair;
        self.cgfunc = cgfunc;
        
        #Initialize the "calculators"
        self.ikfield = None;
        
        if not (not self.__field_obs_names):
            #Construct the field calculator
            self.ikfield = _pyglasstools.LocalCalculator(self.sysdata._getParticleSystem(),self.pair._getPairPotential(),self.cgfunc._getCGFunc())
            
            #Then, add observables
            self.field_obs = observables.initialize_field(self.__field_obs_names,len(self.gridpoints))
            for name in self.field_obs:
                self.ikfield.addObservable(self.field_obs[name]._getObservable())
        
        #Once initialization is done, the logger adds itself to the field list of available loggers 
        loggers.append(self)

        #Finally, we create a list of trajectory files
        #Initialize file to save:
        self.file = {}
        for name in self.__field_obs_names: 
            self.file[name] = MPILogFile(comm, self.fileprefix+"_"+name, MPI.MODE_WRONLY | MPI.MODE_CREATE | MPI.MODE_APPEND)
    def run(self): 
        if not (not self.__field_obs_names):
            self.ikfield.computelocal(self.gridpoints)
    
    def save(self,frame_num):
        for name in self.__field_obs_names:
            val = 0
            signal = False
            if "f_virialstress" in name:
                val = np.array(self.field_obs['f_virialstress'].getVal())[:,int(name[-2]),int(name[-1])]
            else:
                continue

            if rank == 0:
                self.file[name].write("Frame {:d} \n".format(frame_num))
                self.file[name].write("#{:d} \n".format(self.gridsize))
                self.file[name].write("{: <15} {: <15} {: <15} {: <15} \n".format("Coord_x","Coord_y","Coord_z",name))
                for i in range(len(self.gridpoints)):
                    self.file[name].write("{:<15.6e} {:<15.6e} {:<15.6e} {:<15.6e} \n".format(  self.gridpoints[i][0],
                                                                                                self.gridpoints[i][1], 
                                                                                                self.gridpoints[i][2],
                                                                                                val[i]))
                signal = True
                if size > 1:
                    comm.send(signal, dest = rank+1)
            else:
                
                while signal == False:
                    signal = comm.recv(source = rank-1)
                for i in range(len(self.gridpoints)):
                    self.file[name].write("{:<15.6e} {:<15.6e} {:<15.6e} {:<15.6e} \n".format(  self.gridpoints[i][0],
                                                                                                self.gridpoints[i][1], 
                                                                                                self.gridpoints[i][2],
                                                                                                val[i]))
                if rank < size-1:
                    comm.send(signal, dest = rank+1)
    def update(self,frame_num):
        self.sysdata.update(frame_num);
        #Update the field calculator
        self.ikfield = _pyglasstools.LocalCalculator(self.sysdata._getParticleSystem(),self.pair._getPairPotential(),self.cgfunc._getCGFunc())
        #Then, add observables again
        self.field_obs = observables.initialize_field(self.__field_obs_names,len(self.gridpoints))
        for name in self.field_obs:
            self.ikfield.addObservable(self.field_obs[name]._getObservable())

def analyze(frame_list):
    for frame_num in frame_list:
        for logger in loggers:
            logger.update(frame_num);
            logger.run();
            logger.save(frame_num);
