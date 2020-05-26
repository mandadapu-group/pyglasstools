from pyglasstools.nonaffine import _nonaffine
from pyglasstools import _pyglasstools
from pyglasstools import comm, rank, size
import numpy as np

nonaffine_loggers = []

def analyze(frame_list,mode):
    for frame_num in frame_list:
        for logger in nonaffine_loggers:
            logger.update(frame_num);
            logger.run(mode);
            logger.save(frame_num);
        comm.barrier()

#will automatically create the required calculators . . .
class nonaffine_logger(object):
    global loggers 

    def __init__(self, filename=None, names = None, sysdata = None, pair = None, solver = None):
        
        #Save filename
        self.filename = filename
        #Check if we need to track number of converged eigenpairs
        self.getnconv = 'nconv' in names
        #Next, parse the list of names based on what type of obsercables they are
        self.__global_obs_names = [s for s in names if "nonaffine" in s or "eigenvalue" in s];
        self.global_obs = [] 

        self.__eigenvector_obs_names = [s for s in names if "eigenvector" in s ];
        self.eigenvector_obs = [] 
        
        self.__forcedipole_obs_names = [s for s in names if "force_dipole" in s ];
        self.forcedipole_obs = [] 
        #Initialize system data and pair potential of the system
        self.sysdata = sysdata;
        self.pair = pair;
        if(solver == "slepc"):
            self.manager = _nonaffine.PETScManager();
        elif(solver ==  "spectra"):
            self.manager = _nonaffine.HessianManager();
        #Initialize the "hessian" class
        self.hessian = None;
        if not (not self.__global_obs_names) or not (not self.__eigenvector_obs_names) or not( not self.__forcedipole_obs_names):
            self.hessian = hessian(self.sysdata,self.pair,self.manager,comm,solver)
        
        #Once initialization is done, the logger adds itself to the global list of available loggers 
        nonaffine_loggers.append(self)

        #write the initial files
        if not (not self.__global_obs_names) or self.getnconv:
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
        
        if not (not self.__eigenvector_obs_names):
            self.file_eigenvector = {}
            for name in self.__eigenvector_obs_names:
                self.file_eigenvector[name] = _pyglasstools.MPILogFile(comm, name+"_"+solver+".xyz")
        if not (not self.__forcedipole_obs_names):
            self.file_forcedipole = {}
            for name in self.__forcedipole_obs_names:
                self.file_forcedipole[name] = _pyglasstools.MPILogFile(comm, name+"_"+solver+".xyz")

    def run(self,mode="manual"):
        if (mode == "manual"):
            self.hessian.eigs()
            if any("force_dipole" in ext for ext in self.__forcedipole_obs_names):
                self.hessian.solve_forcedipole(1e-2)
        elif (mode=="all"):
            self.hessian.alleigs()
            if any("nonaffine" in ext for ext in self.__global_obs_names):
                self.hessian.compute_nonaffinetensor()
            elif any("force_dipole" in ext for ext in self.__forcedipole_obs_names):
                self.hessian.solve_forcedipole(1e-2)
        elif (mode=="forcedipole"):
            self.hessian.solve_forcedipole(1e-2)
        else:
            raise NameError;
        
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

        if not (not self.__eigenvector_obs_names):
            for name in self.__eigenvector_obs_names:
                if "eigenvector" in name:
                    if rank == 0:
                        self.file_eigenvector[name].write_shared("{:d} \n".format(len(self.sysdata.traj[frame_num].particles.position)))
                        self.file_eigenvector[name].write_shared("#{} {} {} {} {} :".format("Type_ID", "Coord_x","Coord_y","Coord_z",name+"_x",name+"_y"))
                        self.file_eigenvector[name].write_shared("Frame {:d}  \n".format(frame_num))
                    comm.barrier()
                    self.hessian.save_eigenvector(int(name.replace('eigenvector_','')), self.file_eigenvector[name])
                else:
                    continue
        if not (not self.__forcedipole_obs_names):
            for name in self.__forcedipole_obs_names:
                if "force_dipole" in name:
                    if rank == 0:
                        self.file_forcedipole[name].write_shared("{:d} \n".format(len(self.sysdata.traj[frame_num].particles.position)))
                        self.file_forcedipole[name].write_shared("#{} {} {} {} {} :".format("Type_ID", "Coord_x","Coord_y","Coord_z",name+"_x",name+"_y"))
                        self.file_forcedipole[name].write_shared("Frame {:d}  \n".format(frame_num))
                    comm.barrier()
                    self.hessian.save_forcedipole(self.file_forcedipole[name])
                else:
                    continue

    def update(self,frame_num):
        self.sysdata.update(frame_num);
        self.hessian.update(self.sysdata);

class hessian(object):
    def __init__(self, sysdata,potential,manager,comm,solver):
        if (solver == "slepc"):
            self.H = _nonaffine.SLEPcHessian(sysdata._getParticleSystem(),potential._getPairPotential(),manager,comm)
        elif (solver == "spectra"):
            self.H = _nonaffine.SpectraHessian(sysdata._getParticleSystem(),potential._getPairPotential(),manager,comm)
    
    #Redefine attributes so that it directly access Hessian C++ class 
    def __getattr__(self,attr):
            orig_attr = self.H.__getattribute__(attr)
            if callable(orig_attr):
                def hooked(*args, **kwargs):
                    self.pre()
                    result = orig_attr(*args, **kwargs)
                    # prevent H from becoming unwrapped
                    if result == self.H:
                        return self
                    self.post()
                    return result
                return hooked
            else:
                return orig_attr
    def solve_forcedipole(self,threshold):
        return self.H.solveForceDipoleProblem(threshold)
    def save_forcedipole(self,logfile):
        return self.H.saveForceDipoleProblem(logfile)
    def get_eigenvector(self,index):
        return self.H.getEigenvector(index)
    def save_eigenvector(self,index,logfile):
        return self.H.saveEigenvector(index,logfile)
    def get_eigenvalue(self,index):
        return self.H.getEigenvalue(index)
    def save_eigenvalue(self,index,logfile):
        return self.H.saveEigenvalue(index,logfile)
    def get_range(self,index):
        return self.H.getRange(index)
    def update(self,newsysdata):
        self.H.setSystemData(newsysdata._getParticleSystem())
        self.H.buildHessianandMisForce()
    def eigs(self):
        self.H.getEigenPairs()
    def alleigs(self):
        self.H.getAllEigenPairs()
    ## Building nonaffine elasticity tensor
    def compute_nonaffinetensor(self):
        self.H.calculateNonAffineTensor()
    def get_nonaffinetensor(self,i,j, logfile):
        return self.H.saveNonAffineTensor(i,j,logfile)
    def save_nonaffinetensor(self,i,j,logfile):
        return self.H.saveNonAffineTensor(i,j,logfile)
