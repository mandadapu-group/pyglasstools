from pyglasstools.nonaffine import _nonaffine
from pyglasstools import _pyglasstools, comm, rank, size, loggers_list, solvers_list
import numpy as np

def initialize_field(names,dim):
    list_obs = {}
    if any("forcedipole" in s for s in names):
        if dim == 2:
            list_obs['forcedipole'] = _nonaffine.PETScVectorField2D("forcedipole", False,comm)
        elif dim == 3:
            list_obs['forcedipole'] = _nonaffine.PETScVectorField3D("forcedipole", False,comm)
    return list_obs

def initialize_global(names,dim):
    list_obs = {}
    if any("forcedipole_center" in s for s in names):
        if dim == 2:
            list_obs['forcedipole_center'] = _nonaffine.GlobalVector2D("forcedipole_center", "VECTOR", False,comm)
        elif dim == 3:
            list_obs['forcedipole_center'] = _nonaffine.GlobalVector3D("forcedipole_center", "VECTOR", False,comm)
    if any("forcedipole_pj" in s for s in names):
        if dim == 2:
            list_obs['forcedipole_pj'] = _nonaffine.GlobalVector2D("forcedipole_pj", "VECTOR", False,comm)
        elif dim == 3:
            list_obs['forcedipole_pj'] = _nonaffine.GlobalVector3D("forcedipole_pj", "VECTOR", False,comm)
    if any("forcedipole_pi" in s for s in names):
        if dim == 2:
            list_obs['forcedipole_pi'] = _nonaffine.GlobalVector2D("forcedipole_pi", "VECTOR", False,comm)
        elif dim == 3:
            list_obs['forcedipole_pi'] = _nonaffine.GlobalVector3D("forcedipole_pi", "VECTOR", False,comm)
    return list_obs

class logfile(object):
    global loggers_list

    def __init__(self, filename=None, names = None, solver = None, mode="new"):
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
        self.global_obs = initialize_global(self.__obs_names, self.solver.sysdata.simbox.dim)
        self.solver.add_observables(self.global_obs)
        
        Dim = self.solver.sysdata.simbox.dim
        #Create column headers
        if rank == 0 and mode =="new":
            self.file.write_shared("{} ".format("Frame"))
            if not (not self.__obs_names):
                for name in self.__obs_names:
                    if "forcedipole_center" in name or "forcedipole_pi" in name or "forcedipole_pj" in name:
                        for i in range(Dim):
                            self.file.write_shared("{}_{} ".format(name,i))
            self.file.write_shared("\n")
         
    def save(self,frame_num):
        if rank == 0:
            self.file.write_shared("{} ".format(frame_num))
            for name in self.__obs_names:
                Dim = self.solver.sysdata.simbox.dim
                if "forcedipole_center" in name:
                    for i in range(Dim):
                        self.global_obs['forcedipole_center'].save(self.file, i)
                        self.file.write_shared(" ")
                if "forcedipole_pi" in name:
                    for i in range(Dim):
                        self.global_obs['forcedipole_pi'].save(self.file, i)
                        self.file.write_shared(" ")
                if "forcedipole_pj" in name:
                    for i in range(Dim):
                        self.global_obs['forcedipole_pj'].save(self.file, i)
                        self.file.write_shared(" ")
            self.file.write_shared("\n")
    
    def get_forcedipole_center(self):
        vector = np.array(self.global_obs['forcedipole_center'].getValue(),dtype=np.float64)
        if len(vector) < 3:
            vector = np.append(vector,0)
        return vector
class hessian(object):
    global solvers_list
    
    def __init__(self, sysdata,potential, package, checkmineigval=False):

        self.checkmineigval = checkmineigval
        #Initialize system data and pair potential of the system
        self.sysdata = sysdata;
        self.potential = potential;
        if (package == "slepc"):
            self.manager = _nonaffine.PETScManager();
        self.hessian = _nonaffine.SLEPcHessian(sysdata._getParticleSystem(),potential._getPairPotential(),self.manager,comm)
        #elif (package == "spectra"):
        #    self.manager = _nonaffine.HessianManager();
        #    self.hessian = _nonaffine.SpectraHessian(sysdata._getParticleSystem(),potential._getPairPotential(),self.manager,comm)
        self.package = package
        solvers_list.append(self)
    
    def add_observables(self, observables):
        for name in observables:
            if ("forcedipole" in name):
                if "forcedipole_center" in name or "forcedipole_pi" in name or "forcedipole_pj" in name:
                    self.hessian.addGlobalProperty(observables[name])
                else:
                    self.hessian.addVectorField(observables[name])
   
    def runn(self,mode):
        if (mode == "manual"):
            self.eigs()
            if any("force_dipole" in ext for ext in self.__forcedipole_obs_names):
                self.solve_forcedipole()
        elif (mode=="all"):
            self.alleigs()
            if any("nonaffine" in ext for ext in self.__global_obs_names):
                self.compute_nonaffinetensor()
            elif any("force_dipole" in ext for ext in self.__forcedipole_obs_names):
                self.solve_forcedipole()
        elif (mode=="forcedipole"):
            self.solve_forcedipole()
        else:
            raise NameError;
    def run(self,mode="manual"):
        if self.checkmineigval is True:
            self.hessian.checkSmallestEigenvalue_forMumps()
            if self.hessian.mineigval >= 0:
                self.runn(mode)
            else:
                if rank == 0:
                    print("[WARNING] Skip! Negative Eigenvalue detected")
        else:
            self.runn(mode)
    def update(self,frame_num):
        self.sysdata.update(frame_num);
        self.hessian.setSystemData(self.sysdata._getParticleSystem())
        self.hessian.buildHessianandMisForce()
    def set_forcedipole_minmax(forcemin = None, forcemax = None): 
        if forcemin is not None:
           self.manager.fd_random_min = forcemin;     
        if forcemax is not None:
           self.manager.fd_random_max = forcemax;     
    ##Force dipole calculations
    def solve_forcedipole(self):
        if (self.package == "slepc"):
            if (self.manager.fd_mode == "minimum"):
                self.hessian.solveForceDipoleProblem_Minimum()
            elif (self.manager.fd_mode == "random"):
                self.hessian.solveForceDipoleProblem_Random()
        
    #def save_forcedipole(self,logfile):
    #    return self.hessian.saveForceDipoleProblem(logfile)
    #def get_range(self,index):
    #    return self.hessian.getRange(index)
    
    ## computing eigenvectors
    #def eigs(self):
    #    self.hessian.getEigenPairs()
    #def alleigs(self):
    #    self.hessian.getAllEigenPairs()
    #def get_eigenvector(self,index):
    #    return self.hessian.getEigenvector(index)
    #def save_eigenvector(self,index,logfile):
    #    return self.hessian.saveEigenvector(index,logfile)
    #def get_eigenvalue(self,index):
    #    return self.hessian.getEigenvalue(index)
    #def save_eigenvalue(self,index,logfile):
    #    return self.hessian.saveEigenvalue(index,logfile)
    
    ## Building nonaffine elasticity tensor
    #def compute_nonaffinetensor(self):
    #    self.hessian.calculateNonAffineTensor()
    #def get_nonaffinetensor(self,i,j, logfile):
    #    return self.hessian.saveNonAffineTensor(i,j,logfile)
    #def save_nonaffinetensor(self,i,j,logfile):
    #    return self.hessian.saveNonAffineTensor(i,j,logfile)

class fieldlogger(object):
    global loggers_list

    def __init__(self, keyword="dump", names = None, solver = None):
        #Save filename
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
        self.field_obs = initialize_field(self.__obs_names, self.solver.sysdata.simbox.dim)
        self.solver.add_observables(self.field_obs)
        self.file = _pyglasstools.MPILogFile(comm, "{}".format(keyword)+"_"+".xyz")
        
        #Now, we specify the size of particles assigned to a process
        #ends = [sum(counts[:p+1]) for p in range(nprocs)]

    def save_perrank(self,frame_num):
        Dim = self.solver.sysdata.simbox.dim

        ave, res = divmod(len(self.solver.sysdata.traj[frame_num].particles.position), size)
        counts = [ave + 1 if p < res else ave for p in range(size)]
        starts = sum(counts[:rank])
        ends = sum(counts[:rank+1])
        
        for i in range(starts,ends):
            outline = "{} ".format(1)
            outline += "{} ".format(self.solver.sysdata.traj[frame_num].particles.position[i,0])
            outline += "{} ".format(self.solver.sysdata.traj[frame_num].particles.position[i,1])
            if Dim == 3:
                outline += "{} ".format(self.solver.sysdata.traj[frame_num].particles.position[i,2])
            for name in self.__obs_names:
                flattenindex = 0;
                if "forcedipole":
                    outline += self.field_obs["forcedipole"].gettostring(i)
            outline += "\n"
            self.file.write_shared(outline);
    
    def get_vectorfield(self, name,frame_num):
        vector = []
        ave, res = divmod(len(self.solver.sysdata.traj[frame_num].particles.position), size)
        counts = [ave + 1 if p < res else ave for p in range(size)]
        starts = sum(counts[:rank])
        ends = sum(counts[:rank+1])
        for i in range(starts,ends):
            vector.append(np.reshape(self.field_obs[name].getVectorValue(i),(3,)))
        #Next, we gather all of these vectors
        newvector = comm.all_gather_v(vector)
        vector = np.reshape(newvector,(len(self.solver.sysdata.traj[frame_num].particles.position),3));
        del newvector;
        return vector#Reshape
        
    def save(self,frame_num):
        if rank == 0:
            self.file.write_shared("{:d} \n".format(len(self.solver.sysdata.traj[frame_num].particles.position)))
            self.file.write_shared("#Frame {:d}  \n".format(frame_num))
        comm.barrier()
        self.save_perrank(frame_num)
