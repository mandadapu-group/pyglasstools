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
    if any("eigenvector" in s for s in names):
        for s in [name for name in names if "eigenvector" in name]:
            if dim == 2:
                list_obs[s] = _nonaffine.PETScVectorField2D(s, False,comm)
            elif dim == 3:
                list_obs[s] = _nonaffine.PETScVectorField3D(s, False,comm)
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
    if any("nonaffinetensor" in s for s in names):
        if dim == 2:
            list_obs['nonaffinetensor'] = _nonaffine.NonAffineTensor2D("nonaffinetensor", "4-TENSOR", False,comm)
    if any("eigenvalue" in s for s in names):
        for s in [name for name in names if "eigenvalue" in name]:
            list_obs[s] = _nonaffine.GlobalScalar(s, "SCALAR", False,comm)
    return list_obs

class fdcalculator(object):
    global solvers_list

    def __init__(self, sysdata, potential):
        self.pyhessian = hessian(sysdata,potential,"petsc");
        self.cppcalculator = _nonaffine.PETScForceDipoleCalculator(self.pyhessian.cpphessian)
        solvers_list.append(self)

    def add_observables(self, observables):
        for name in observables:
            if ("forcedipole" in name):
                if "forcedipole_center" in name or "forcedipole_pi" in name or "forcedipole_pj" in name:
                    self.cppcalculator.addGlobalProperty(observables[name])
                else:
                    self.cppcalculator.addVectorField(observables[name])
    
    def run(self,mode="manual"):
        self.cppcalculator.solveForceDipoleProblem()
    
    def set_forcedipole_minmax(forcemin = None, forcemax = None): 
        if forcemin is not None:
           self.pyhessian.cppmanager.fd_random_min = forcemin;     
        if forcemax is not None:
           self.pyhessian.cppmanager.fd_random_max = forcemax;     
    
    def update(self,frame_num):
        self.pyhessian.update(frame_num)
        self.cppcalculator.setHessian(self.pyhessian.cpphessian)
        
class eigensolver(object):
    global solvers_list

    def __init__(self, sysdata, potential, package = "slepc"):
        if package == "slepc":
            self.pyhessian = hessian(sysdata,potential,"slepc");
            self.cppeigensolver = _nonaffine.SLEPcNMA(self.pyhessian.cpphessian)
        #We need Spectra implementation here as well
        solvers_list.append(self)

    def add_observables(self, observables):
        for name in observables:
            if "nonaffinetensor" in name:
                self.cppeigensolver.addGlobalProperty(observables[name])
            if "eigenvalue" in name:
                self.cppeigensolver.addGlobalProperty(observables[name])
            if "eigenvector" in name:
                self.cppeigensolver.addVectorField(observables[name])
    
    def run(self,mode="manual"):
        self.cppeigensolver.getAllEigenPairs()
    
    def update(self,frame_num):
        self.pyhessian.update(frame_num)
        self.cppeigensolver.setHessian(self.pyhessian.cpphessian)

class hessian(object):
    def __init__(self, sysdata,potential, package):
        #Initialize system data and pair potential of the system
        self.sysdata = sysdata;
        self.potential = potential;
        self.cppmanager = _nonaffine.PETScManager();
        if (package == "petsc"):
            self.cpphessian = _nonaffine.PETScHessian2D(sysdata.particledata,potential._getPairPotential(),self.cppmanager,comm)
        elif (package == "slepc"):
            self.cpphessian = _nonaffine.SLEPcHessian2D(sysdata.particledata,potential._getPairPotential(),self.cppmanager,comm)
        self.frame_num = self.sysdata.frame_num
        #elif (package == "spectra"):
        #    self.manager = _nonaffine.HessianManager();
        #    self.hessian = _nonaffine.SpectraHessian(sysdata._getParticleSystem(),potential._getPairPotential(),self.manager,comm)
        self.package = package
   
    def update(self,frame_num):
        self.sysdata.update(frame_num);
        self.cpphessian.setSystemData(self.sysdata.particledata)
        if self.frame_num != frame_num:
            self.cpphessian.destroyPETScObjects()
            self.cpphessian.assemblePETScObjects()
            self.frame_num = frame_num


########

class logfile(object):
    global loggers_list

    def __init__(self, filename=None, names = None, solver = None, sysdata = None, mode="new"):
        #Save filename
        self.filename = filename
        #Next, parse the list of names based on what type of obsercables they are
        self.__obs_names = names
        #[s for s in names if "" in s or "eigenvalue" in s];
        self.sysdata = sysdata
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
        self.global_obs = initialize_global(self.__obs_names, self.sysdata.simbox.dim)
        self.solver.add_observables(self.global_obs)
        
        Dim = self.sysdata.simbox.dim
        #Create column headers
        if rank == 0 and mode =="new":
            self.file.write_shared("{} ".format("Frame"))
            if not (not self.__obs_names):
                for name in self.__obs_names:
                    if "forcedipole_center" in name or "forcedipole_pi" in name or "forcedipole_pj" in name:
                        for i in range(Dim):
                            self.file.write_shared("{}_{} ".format(name,i))
                    else:
                            self.file.write_shared("{} ".format(name))
            self.file.write_shared("\n")
         
    def save(self,frame_num):
        if rank == 0:
            self.file.write_shared("{} ".format(frame_num))
            for name in self.__obs_names:
                Dim = self.sysdata.simbox.dim
                if "forcedipole" in name:
                    for i in range(Dim):
                        self.global_obs[name].save(self.file, i)
                        self.file.write_shared(" ")
                if "eigenvalue" in name:
                    self.global_obs[name].save(self.file)
                if "nonaffinetensor" in name:
                    i = int(name[-4]) 
                    j = int(name[-3])
                    k = int(name[-2]) 
                    l = int(name[-1])
                    self.global_obs['nonaffinetensor'].save(self.file,l+Dim*(k+Dim*(j+Dim*i)))
                self.file.write_shared(" ")
            self.file.write_shared("\n")
    
    def get_forcedipole_center(self):
        vector = np.array(self.global_obs['forcedipole_center'].getValue(),dtype=np.float64)
        if len(vector) < 3:
            vector = np.append(vector,0)
        return vector

class fieldlogger(object):
    global loggers_list

    def __init__(self, keyword="dump", names = None, solver = None, sysdata = None):
        #Save filename
        self.keyword = keyword
        #Next, parse the list of names based on what type of obsercables they are
        self.__obs_names = names
        #[s for s in names if "" in s or "eigenvalue" in s];
        self.sysdata = sysdata;
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
        self.field_obs = initialize_field(self.__obs_names, self.sysdata.simbox.dim)
        self.solver.add_observables(self.field_obs)
        self.file = _pyglasstools.MPILogFile(comm, "{}".format(keyword)+"_"+".xyz")
        
        #Now, we specify the size of particles assigned to a process
        #ends = [sum(counts[:p+1]) for p in range(nprocs)]

    def save_perrank(self,frame_num):
        Dim = self.sysdata.simbox.dim

        ave, res = divmod(len(self.sysdata.traj[frame_num].particles.position), size)
        counts = [ave + 1 if p < res else ave for p in range(size)]
        starts = sum(counts[:rank])
        ends = sum(counts[:rank+1])
        
        for i in range(starts,ends):
            outline = "{} ".format(1)
            outline += "{} ".format(self.sysdata.traj[frame_num].particles.position[i,0])
            outline += "{} ".format(self.sysdata.traj[frame_num].particles.position[i,1])
            if Dim == 3:
                outline += "{} ".format(self.sysdata.traj[frame_num].particles.position[i,2])
            for name in self.__obs_names:
                flattenindex = 0;
                if "forcedipole" in name:
                    outline += self.field_obs["forcedipole"].gettostring(i)
                if "eigenvector" in name:
                    outline += self.field_obs[name].gettostring(i)
            outline += "\n"
            self.file.write_shared(outline);
    
    def get_vectorfield(self, name,frame_num):
        vector = []
        ave, res = divmod(len(self.sysdata.traj[frame_num].particles.position), size)
        counts = [ave + 1 if p < res else ave for p in range(size)]
        starts = sum(counts[:rank])
        ends = sum(counts[:rank+1])
        for i in range(starts,ends):
            vector.append(np.reshape(self.field_obs[name].getVectorValue(i),(3,)))
        #Next, we gather all of these vectors
        newvector = comm.all_gather_v(vector)
        vector = [item for sublist in newvector for item in sublist]
        #vector = np.reshape(newvector,(len(self.sysdata.traj[frame_num].particles.position),3));
        del newvector;
        return vector#Reshape
        
    def save(self,frame_num):
        if rank == 0:
            self.file.write_shared("{:d} \n".format(len(self.sysdata.traj[frame_num].particles.position)))
            self.file.write_shared("#Frame {:d}  \n".format(frame_num))
        comm.barrier()
        self.save_perrank(frame_num)
