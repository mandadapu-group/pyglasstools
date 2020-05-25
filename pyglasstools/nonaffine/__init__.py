from pyglasstools.nonaffine import _nonaffine
from pyglasstools import _pyglasstools
from pyglasstools import MPILogFile
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
size = nprocs
nonaffine_loggers = []

def analyze(frame_list,mode):
    for frame_num in frame_list:
        for logger in nonaffine_loggers:
            logger.update(frame_num);
            logger.run(mode);
            logger.save(frame_num);
        comm.Barrier()

#will automatically create the required calculators . . .
class slepc_logger(object):
    global loggers 

    def __init__(self, filename=None, names = None, sysdata = None, pair = None):
        
        #Save filename
        self.filename = filename
        
        #Check if we need to track number of converged eigenpairs
        self.getnconv = 'nconv' in names

        #Next, parse the list of names based on what type of obsercables they are
        self.__global_obs_names = [s for s in names if "nonaffine" in s or "eigenvalue" in s];
        self.global_obs = [] 

        self.__eigenvector_obs_names = [s for s in names if "eigenvector" in s ];
        self.eigenvector_obs = [] 
        
        #Initialize system data and pair potential of the system
        self.sysdata = sysdata;
        self.pair = pair;
        self.manager = _nonaffine.PETScManager();
        
        #Initialize the "hessian" class
        self.hessian = None;
        self.communicator = _pyglasstools.Communicator(); 
        if not (not self.__global_obs_names) or not (not self.__eigenvector_obs_names):
            #Construct the hessian
            self.hessian = slepc_hessian(self.sysdata,self.pair,self.manager)
        
        #Once initialization is done, the logger adds itself to the global list of available loggers 
        nonaffine_loggers.append(self)
        

        #write the initial files
        if not (not self.__global_obs_names) or self.getnconv:
            #Initialize file to save:
            self.file = MPILogFile(comm, self.filename, MPI.MODE_WRONLY | MPI.MODE_CREATE | MPI.MODE_APPEND)

            #Create column headers
            if rank == 0:
                self.file.write("{:<10} ".format("Frame"))
                if self.getnconv:
                    self.file.write("{:<10} ".format('nconv'))
                if not (not self.__global_obs_names):
                    for name in self.__global_obs_names:
                        self.file.write("{:<20} ".format(name))
                
                self.file.write("\n")
        
        if not (not self.__eigenvector_obs_names):
            self.file_eigenvector = {}
            for name in self.__eigenvector_obs_names:
                self.file_eigenvector[name] = MPILogFile(comm, name+".xyz", MPI.MODE_WRONLY | MPI.MODE_CREATE | MPI.MODE_APPEND)

    def run(self,mode="manual"):
        if (mode == "manual"):
            self.hessian.eigs()
        elif (mode=="all"):
            self.hessian.alleigs()
        else:
            raise NameError;
        if any("nonaffine" in ext for ext in self.__global_obs_names):
            self.hessian.compute_nonaffine()
        #else:
        #    self.hessian.eigs()

    def save(self,frame_num):
        if not (not self.__global_obs_names):
            if rank == 0:
                self.file.write("{:<20d} ".format(frame_num))
                if self.getnconv:
                    self.file.write("{:<10d} ".format(self.hessian.nconv))

            for name in self.__global_obs_names:
                val = 0
                if "nonaffine" in name:
                    val = self.hessian.nonaffinetensor[int(name[-2]),int(name[-1])]
                elif "eigenvalue" in name:
                    val = self.hessian.get_eigenvalue(int(name.replace('eigenvalue_','')))
                
                if rank == 0:
                    self.file.write("{:<20.12e} ".format(val))
            if rank == 0:
                self.file.write("\n")
        comm.Barrier()

        if not (not self.__eigenvector_obs_names):
            for name in self.__eigenvector_obs_names:
                val = 0
                if "eigenvector" in name:
                    #check the name
                    val = np.array(self.hessian.get_eigenvector(int(name.replace('eigenvector_',''))))
                else:
                    continue
                if rank == 0:
                    self.file_eigenvector[name].write("{:d} \n".format(len(self.sysdata.traj[frame_num].particles.position)))
                    self.file_eigenvector[name].write("#{: <14} {: <15} {: <15} {: <15} :".format("Coord_x","Coord_y","Coord_z",name+"_x",name+"_y",name+"_mag"))
                    self.file_eigenvector[name].write("Frame {:d}  \n".format(frame_num))
                    for i in range(len(self.sysdata.traj[frame_num].particles.position)):
                        self.file_eigenvector[name].write("{:d} {:<15.6e} {:<15.6e} {:<15.6e} {:<15.6e} {:<15.6e} {:<15.6e} \n".format(   self.sysdata.traj[frame_num].particles.typeid[i],
                                                                                                                                self.sysdata.traj[frame_num].particles.position[i,0],
                                                                                                                                self.sysdata.traj[frame_num].particles.position[i,1],
                                                                                                                                self.sysdata.traj[frame_num].particles.position[i,2],
                                                                                                                                val[2*i],val[2*i+1],np.sqrt(val[2*i]**2+val[2*i+1]**2)))
                comm.Barrier();

    def update(self,frame_num):
        self.sysdata.update(frame_num);
        self.hessian.update(self.sysdata);

class spectra_logger(object):
    global loggers 

    def __init__(self, filename=None, names = None, sysdata = None, pair = None, solver= None):
        
        #Save filename
        self.filename = filename
        
        self.getnconv = 'nconv' in names
        #First, parse the list of names based on what type of obsercables they are
        self.__global_obs_names = [s for s in names if "nonaffine" in s or "eigenvalue" in s];
        self.global_obs = [] 

        self.__eigenvector_obs_names = [s for s in names if "eigenvector" in s ];
        self.eigenvector_obs = [] 
        
        #Initialize system data and pair potential of the system
        self.sysdata = sysdata;
        self.pair = pair;
        self.manager = _nonaffine.HessianManager();
        
        #Initialize the "hessian" class
        self.hessian = None;
        self.communicator = _pyglasstools.Communicator(); 
        self.solver = solver
        if not (not self.__global_obs_names) or not (not self.__eigenvector_obs_names):
            #Construct the hessian
            self.hessian = spectra_hessian(self.sysdata,self.pair,self.communicator,self.manager)
        
        #Once initialization is done, the logger adds itself to the global list of available loggers 
        nonaffine_loggers.append(self)

        if not (not self.__global_obs_names) or self.getnconv:
            #Initialize file to save:
            self.file = MPILogFile(comm, self.filename, MPI.MODE_WRONLY | MPI.MODE_CREATE | MPI.MODE_APPEND)

            #Create column headers
            if rank == 0:
                self.file.write("{:<10} ".format("Frame"))
                if self.getnconv:
                    self.file.write("{:<10} ".format('nconv'))
                if not (not self.__global_obs_names):
                    for name in self.__global_obs_names:
                        self.file.write("{:<20} ".format(name))
                
                self.file.write("\n")
        
        if not (not self.__eigenvector_obs_names):
            self.file_eigenvector = {}
            for name in self.__eigenvector_obs_names:
                self.file_eigenvector[name] = MPILogFile(comm, name+".xyz", MPI.MODE_WRONLY | MPI.MODE_CREATE | MPI.MODE_APPEND)

    def run(self,mode="manual"):
        if (mode == "manual"):
            self.hessian.eigs()
        elif (mode=="all"):
            self.hessian.alleigs()
        else:
            raise NameError;
        if any("nonaffine" in ext for ext in self.__global_obs_names):
            self.hessian.compute_nonaffine()
        #else:
        #    self.hessian.eigs()

    def save(self,frame_num):
        if not (not self.__global_obs_names):
            if rank == 0:
                self.file.write("{:<20d} ".format(frame_num))
                if self.getnconv:
                    self.file.write("{:<10d} ".format(self.hessian.nconv))

            for name in self.__global_obs_names:
                if "nonaffine" in name:
                    val = self.hessian.get_nonaffinetensor(int(name[-2]),int(name[-1]))
                elif "eigenvalue" in name:
                    val = self.hessian.get_eigenvalue(int(name.replace('eigenvalue_','')))
                if val is not None:
                    self.file.write("{:<20.12e} ".format(val))
            if rank == 0:
                self.file.write("\n")
        comm.Barrier()

        if not (not self.__eigenvector_obs_names):
            for name in self.__eigenvector_obs_names:
                val = 0
                if "eigenvector" in name:
                    #check the name
                    val = np.array(self.hessian.get_eigenvector(int(name.replace('eigenvector_',''))))
                else:
                    continue
                if rank == 0:
                    self.file_eigenvector[name].write("{:d} \n".format(len(self.sysdata.traj[frame_num].particles.position)))
                    self.file_eigenvector[name].write("#{: <14} {: <15} {: <15} {: <15} :".format("Coord_x","Coord_y","Coord_z",name+"_x",name+"_y",name+"_mag"))
                    self.file_eigenvector[name].write("Frame {:d}  \n".format(frame_num))
                    for i in range(len(self.sysdata.traj[frame_num].particles.position)):
                        self.file_eigenvector[name].write("{:d} {:<15.6e} {:<15.6e} {:<15.6e} {:<15.6e} {:<15.6e} {:<15.6e} \n".format(   self.sysdata.traj[frame_num].particles.typeid[i],
                                                                                                                                self.sysdata.traj[frame_num].particles.position[i,0],
                                                                                                                                self.sysdata.traj[frame_num].particles.position[i,1],
                                                                                                                                self.sysdata.traj[frame_num].particles.position[i,2],
                                                                                                                                val[2*i],val[2*i+1],np.sqrt(val[2*i]**2+val[2*i+1]**2)))
                comm.Barrier();

    def update(self,frame_num):
        self.sysdata.update(frame_num);
        self.hessian.update(self.sysdata);

class slepc_hessian(object):
    def __init__(self, sysdata,potential,manager):
        self.H = _nonaffine.SLEPcHessian(sysdata._getParticleSystem(),potential._getPairPotential(),manager)
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
    def get_eigenvector(self,index):
        return self.H.getEigenvector(index)
    
    def get_eigenvalue(self,index):
        return self.H.getEigenvalue(index)

    def get_range(self,index):
        return self.H.getRange(index)
    def update(self,newsysdata):
        self.H.setSystemData(newsysdata._getParticleSystem())
        self.H.buildHessianandMisForce()
    def eigs(self):
        self.H.getEigenPairs()
    def alleigs(self):
        self.H.getAllEigenPairs_Mumps()
    ## Building nonaffine elasticity tensor
    def compute_nonaffine(self):
        self.H.calculateNonAffine()
            #Merge the computed pseudoinverses from each process
            #And store to Root Process
            #self.assembled_nonaffine = comm.reduce(self.H.nonaffinetensor,MPI.SUM,root=0)

class spectra_hessian(object):
    def __init__(self, sysdata,potential,comm,manager):
        self.H = _nonaffine.Hessian(sysdata._getParticleSystem(),potential._getPairPotential(), manager, comm)
     
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
    
    ## Our Implementation of eigendecomposition and pseudoinverse
    def eigs(self, selrule = 'LM', nev = 1, ncv = 5, maxiter=1000, tol=1e-10):
        self.H.getEigenDecomposition(selrule,nev,ncv,maxiter);
    ## Shortcut Function to Compute All Eigenpairs
    def alleigs(self,dim=2,maxiter=10000, tol=1e-8):
        self.H.getEigenDecomposition('LM',self.H.hessian.shape[0]-dim,self.H.hessian.shape[0],maxiter);
    def get_eigenvector(self,index):
        return self.H.getEigenvector(index)
    def get_nonaffinetensor(self,i,j):
        val =  self.H.getNonAffinetensor(i,j)
        if val[0]:
            return val[1]
        else:
            return None
    def get_eigenvalue(self,index):
        val =  self.H.getEigenvalue(index)
        if val[0]:
            return val[1]
        else:
            return None
    def update(self,newsysdata):
        self.H.setSystemData(newsysdata._getParticleSystem())
        self.H.buildHessianandMisForce()
    
    ## Building pseudoinverse
    def build_pinv(self,tol):
        self.H.buildPseudoInverse(tol)
    
    ## Building nonaffine elasticity tensor
    def compute_nonaffine(self):
        self.H.calculateNonAffine()
    
    def _getObservable(self):
        return self.H
