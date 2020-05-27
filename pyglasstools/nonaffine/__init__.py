from pyglasstools.nonaffine import _nonaffine
from pyglasstools import _pyglasstools
from pyglasstools import comm, rank, size
import numpy as np

loggers_list = []
solvers_list = []
def analyze(frame_list,mode):
    for frame_num in frame_list:
        for solver in solvers_list:
            solver.update(frame_num);
            solver.run(mode);
        comm.barrier()
        for logger in loggers_list:
            logger.save(frame_num);

class logfile(object):
    global loggers_list

    def __init__(self, filename=None, names = None, solver = None):
        #Save filename
        self.filename = filename
        
        #Next, parse the list of names based on what type of obsercables they are
        self.__global_obs_names = [s for s in names if "nonaffine" in s or "eigenvalue" in s];
        self.global_obs = [] 

        #Initialize the "hessian" class
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
        if not (not self.__global_obs_names):
            #Initialize file to save:
            self.file = _pyglasstools.MPILogFile(comm, self.filename)

            #Create column headers
            if rank == 0:
                self.file.write_shared("{} ".format("Frame"))
                if self.getnconv:
                    self.file.write_shared("{} ".format('nconv'))
                if not (not self.__global_obs_names):
                    for name in self.__global_obs_names:
                        self.file.write_shared("{} ".format(name))
                self.file.write_shared("\n")
    
    def save(self,frame_num):
        if not (not self.__global_obs_names):
            if rank == 0:
                self.file.write_shared("{:<20d} ".format(frame_num))
                if self.getnconv:
                    self.file.write_shared("{:<10d} ".format(self.hessian.nconv))
            
            for name in self.__global_obs_names:
                if "nonaffine" in name:
                    self.hessian.save_nonaffinetensor(int(name[-2]),int(name[-1]), self.file)
                elif "eigenvalue" in name:
                    self.hessian.save_eigenvalue(int(name.replace('eigenvalue_','')),self.file)
                if rank == 0:
                    self.file.write_shared("\n")

#will automatically create the required calculators . . .
class fieldlog(object):
    global loggers_list 

    def __init__(self, keyword=None, names = None, solver = None):
        #Save filename
        self.keyword = keyword
        if solver is None:
            if rank == 0:
                print("[ERROR] User must specify solver object!")
            comm.abort(1)
        else:
            self.solver = solver;

        if names is None:
            if rank == 0:
                print("[ERROR] User must specify list of observables to save!")
            comm.abort(1)
        else:
            self.__obs_names = [s for s in names if "eigenvector" or "force_dipole" in s ];
            if not self.__obs_names:
                print("[ERROR] Observables are not recognized. Only eigenvectors and force_dipole are available!")
                comm.abort(1)
        #Once initialization is done, the logger adds itself to the global list of available loggers 
        loggers_list.append(self)
            
        self.file_list = {}
        for name in self.__obs_names:
            self.file_list[name] = _pyglasstools.MPILogFile(comm, "{}".format(keyword)+name+"_"+".xyz")
        
    def save(self,frame_num):
        for name in self.__obs_names:
            if rank == 0:
                self.file_list[name].write_shared("{:d} \n".format(len(self.solver.sysdata.traj[frame_num].particles.position)))
                self.file_list[name].write_shared("#Frame {:d}  \n".format(frame_num))
            if "eigenvector" in name:
                self.solver.save_eigenvector(int(name.replace('eigenvector_','')), self.file_list[name])
            elif "force_dipole" in name:
                self.solver.save_forcedipole(self.file_list[name])

class hessian(object):
    global solvers_list
    def __init__(self, sysdata,potential, package):

        #Initialize system data and pair potential of the system
        self.sysdata = sysdata;
        self.potential = potential;
        if (package == "slepc"):
            self.manager = _nonaffine.PETScManager();
            self.hessian = _nonaffine.SLEPcHessian(sysdata._getParticleSystem(),potential._getPairPotential(),self.manager,comm)
        elif (package == "spectra"):
            self.manager = _nonaffine.HessianManager();
            self.hessian = _nonaffine.SpectraHessian(sysdata._getParticleSystem(),potential._getPairPotential(),self.manager,comm)
        self.package = package
        solvers_list.append(self)
    def run(self,mode="manual"):
        if (mode == "manual"):
            self.eigs()
            if any("force_dipole" in ext for ext in self.__forcedipole_obs_names):
                self.solve_forcedipole(1e-4)
        elif (mode=="all"):
            self.alleigs()
            if any("nonaffine" in ext for ext in self.__global_obs_names):
                self.compute_nonaffinetensor()
            elif any("force_dipole" in ext for ext in self.__forcedipole_obs_names):
                self.solve_forcedipole(1e-4)
        elif (mode=="forcedipole"):
            self.solve_forcedipole(1e-4)
        else:
            raise NameError;
    
    def update(self,frame_num):
        self.sysdata.update(frame_num);
        self.hessian.setSystemData(self.sysdata._getParticleSystem())
        self.hessian.buildHessianandMisForce()
   
    ##Force dipole calculations
    def solve_forcedipole(self,threshold):
        if (self.package == "slepc"):
            return self.hessian.solveForceDipoleProblem(threshold)
        else:
            return None
    def save_forcedipole(self,logfile):
        return self.hessian.saveForceDipoleProblem(logfile)
    def get_range(self,index):
        return self.hessian.getRange(index)
    
    ## computing eigenvectors
    def eigs(self):
        self.hessian.getEigenPairs()
    def alleigs(self):
        self.hessian.getAllEigenPairs()
    def get_eigenvector(self,index):
        return self.hessian.getEigenvector(index)
    def save_eigenvector(self,index,logfile):
        return self.hessian.saveEigenvector(index,logfile)
    def get_eigenvalue(self,index):
        return self.hessian.getEigenvalue(index)
    def save_eigenvalue(self,index,logfile):
        return self.hessian.saveEigenvalue(index,logfile)
    
    ## Building nonaffine elasticity tensor
    def compute_nonaffinetensor(self):
        self.hessian.calculateNonAffineTensor()
    def get_nonaffinetensor(self,i,j, logfile):
        return self.hessian.saveNonAffineTensor(i,j,logfile)
    def save_nonaffinetensor(self,i,j,logfile):
        return self.hessian.saveNonAffineTensor(i,j,logfile)
