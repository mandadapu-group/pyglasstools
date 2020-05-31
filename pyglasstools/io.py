from pyglasstools import _pyglasstools
from pyglasstools import thermo
from pyglasstools import irvingkirkwood
from pyglasstools import comm, rank, size, loggers_list
import numpy as np

class logfile(object):
    global loggers_list

    def __init__(self, filename=None, names = None, solver = None):
        #Save filename
        self.filename = filename
        #Next, parse the list of names based on what type of obsercables they are
        self.__obs_names = names
        #[s for s in names if "" in s or "eigenvalue" in s];

        #Initialize a thermoproperty class
        if solver is None:
            if rank == 0:
                print("[ERROR] User must specify solver object!")
            comm.abort(1)
        elif names is None:
            if rank == 0:
                print("[ERROR] User must specify list of observables to save!")
            comm.abort(1)
        elif filename is None:
            if rank == 0:
                print("[ERROR] User must specify a filename for the logfile!")
            comm.abort(1)
        
        self.solver = solver;
        
        #Once initialization is done, the logger adds itself to the global list of available loggers 
        loggers_list.append(self)

        #write the initial files
        #Initialize file to save:
        self.file = _pyglasstools.MPILogFile(comm, self.filename)
        
        #Then, add observables
        self.global_obs = thermo.initialize_global(self.__obs_names, self.solver.sysdata.simbox.dim)
       
        self.solver.add_observables(self.global_obs)
        #Create column headers
        if rank == 0:
            self.file.write_shared("{} ".format("Frame"))
            if not (not self.__obs_names):
                for name in self.__obs_names:
                    self.file.write_shared("{} ".format(name))
            self.file.write_shared("\n")
         
    def save(self,frame_num):
        if rank == 0:
            self.file.write_shared("{} ".format(frame_num))
        for name in self.__obs_names:
            Dim = self.solver.sysdata.simbox.dim
            if "virialstress" in name:
                i = int(name[-2]) 
                j = int(name[-1])
                self.global_obs['virialstress'].save(self.file, i+Dim*j)
                self.file.write_shared(" ")
            elif "kineticstress" in name:
                i = int(name[-2]) 
                j = int(name[-1])
                self.global_obs['kineticstress'].save(self.file, i+Dim*j)
                self.file.write_shared(" ")
            elif "borntensor" in name:
                i = int(name[-4]) 
                j = int(name[-3])
                k = int(name[-2]) 
                l = int(name[-1])
                #self.global_obs['borntensor'].save(self.file,i+Dim*(j+Dim*(k+Dim*l)))
                self.global_obs['borntensor'].save(self.file,l+Dim*(k+Dim*(j+Dim*i)))
                self.file.write_shared(" ")
        self.file.write_shared("\n")

class fieldlogger(object):
    global loggers_list

    def __init__(self, keyword=None, names = None, solver = None):
        #Save filename
        if (keyword == None):
            self.keyword = ""
        else:
            self.keyword = keyword
        #Next, parse the list of names based on what type of obsercables they are
        self.__obs_names = names
        #[s for s in names if "" in s or "eigenvalue" in s];

        #Initialize a thermoproperty class
        if solver is None:
            if rank == 0:
                print("[ERROR] User must specify solver object!")
            comm.abort(1)
        elif names is None:
            if rank == 0:
                print("[ERROR] User must specify list of observables to save!")
            comm.abort(1)
        
        self.solver = solver;
        
        #Once initialization is done, the logger adds itself to the global list of available loggers 
        loggers_list.append(self)

        #write the initial files
        #Then, add observables
        self.field_obs = irvingkirkwood.initialize_field(self.__obs_names, self.solver.sysdata.simbox.dim,self.solver.gridpoints)
        self.solver.add_observables(self.field_obs)
        self.file = _pyglasstools.MPILogFile(comm, "{}".format(keyword)+"_"+".xyz")
         
    def save(self,frame_num):
        Dim = self.solver.sysdata.simbox.dim
        if rank == 0:
            self.file.write_shared("{:d} \n".format(len(self.solver.gridpoints)))
            self.file.write_shared("#Frame {:d}  \n".format(frame_num))
        for index, gridpos in enumerate(self.solver.gridpoints):
            self.file.write_shared("{} ".format(self.solver.gridpoints[index][0]))
            self.file.write_shared("{} ".format(self.solver.gridpoints[index][1]))
            if Dim == 3:
                self.file.write_shared("{} ".format(self.solver.gridpoints[index][2]))
            for name in self.__obs_names:
                flattenindex = 0;
                if "stress" in name:
                    i = int(name[-2]) 
                    j = int(name[-1])
                    flattenindex = i+Dim*j
                self.field_obs[name[:-3]].save(self.file, i+Dim*j, index)
            self.file.write_shared("\n");
