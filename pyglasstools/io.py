from pyglasstools import _pyglasstools
from pyglasstools import thermo
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
            if "virialstress" in name:
                self.global_obs['virialstress'].save(self.file, [int(name[-2]),int(name[-1])])
                self.file.write_shared(" ")
            elif "kineticstress" in name:
                self.global_obs['virialstress'].save(self.file, [int(name[-2]),int(name[-1])])
                self.file.write_shared(" ")
            elif "borntensor" in name:
                self.global_obs['virialstress'].save(self.file, [int(name[-2]),int(name[-1])])
                self.file.write_shared(" ")
        self.file.write_shared("\n")
