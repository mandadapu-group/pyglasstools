import pyglasstools
from pyglasstools.elasticity import _elasticity
from . import _elasticitypetsc
from pyglasstools import _pyglasstools, comm, rank, size, loggers_list, solvers_list
import numpy as np
import os

def initialize_field(names,dim):
    list_obs = {}
    if any("displacement" in s for s in names):
        if dim == 2:
            list_obs['displacement'] = _elasticitypetsc.PETScVectorField2D("displacement", False,comm)
        elif dim == 3:
            list_obs['displacement'] = _elasticitypetsc.PETScVectorField3D("displacement", False,comm)
    if any("eigenvector" in s for s in names):
        for s in [name for name in names if "eigenvector" in name]:
            if dim == 2:
                list_obs[s] = _elasticitypetsc.PETScVectorField2D(s, False,comm)
            elif dim == 3:
                list_obs[s] = _elasticitypetsc.PETScVectorField3D(s, False,comm)
    if any("loclandscape" in s for s in names):
        if dim == 2:
            list_obs["loclandscape"] = _elasticitypetsc.PETScVectorField2D("loclandscape", False,comm)
        elif dim == 3:
            list_obs["loclandscape"] = _elasticitypetsc.PETScVectorField3D("loclandscape", False,comm)
    return list_obs

def initialize_global(names,dim):
    list_obs = {}
    if any("nonaffinetensor" in s for s in names):
        if dim == 2:
            list_obs['nonaffinetensor'] = _elasticity.NonAffineTensor2D("nonaffinetensor", "4-TENSOR", False,comm)
        elif dim == 3:
            list_obs['nonaffinetensor'] = _elasticity.NonAffineTensor3D("nonaffinetensor", "4-TENSOR", False,comm)
    if any("eigenvalue" in s for s in names):
        for s in [name for name in names if "eigenvalue" in name]:
            list_obs[s] = _elasticity.GlobalScalar(s, "SCALAR", False, dim, comm)
    if any("eigenrelerror" in s for s in names):
        for s in [name for name in names if "eigenrelerror" in name]:
            list_obs[s] = _elasticity.GlobalScalar(s, "SCALAR", False, dim, comm)
    if any("nconv" in s for s in names):
        list_obs["nconv"] = _elasticity.GlobalScalar("nconv", "SCALAR", False, dim, comm)
    return list_obs

class linrescalculator(object):
    global solvers_list

    def __init__(self):
        self.pyhessian = hessian("petsc");
        self.cppcalculator = _elasticitypetsc.PETScLinearResponse(self.pyhessian.cpphessian)
        solvers_list.append(self)

    def add_observables(self, observables):
        for name in observables:
            if ("displacement" in name):
                self.cppcalculator.addVectorField(observables[name])
    
    def run(self):
        self.cppcalculator.solveLinearResponseProblem()
    
    def update(self,frame_num):
        self.pyhessian.update(frame_num)
        self.cppcalculator.setHessian(self.pyhessian.cpphessian)
        
class eigensolver(object):
    global solvers_list

    def __init__(self, package = "slepc-petsc"):
        self.package = package
        if package == "slepc-petsc" or package == "slepc-mumps":
            self.pyhessian = hessian("slepc");
            self.cppeigensolver = _elasticitypetsc.SLEPcNMA(self.pyhessian.cpphessian)
        solvers_list.append(self)

    def add_observables(self, observables):
        for name in observables:
            if "nonaffinetensor" in name or "nconv" in name or "eigenvalue" in name or "eigenrelerror" in name:
                self.cppeigensolver.addGlobalProperty(observables[name])
            if "eigenvector" in name:
                self.cppeigensolver.addVectorField(observables[name])
            if "loclandscape" in name:
                self.cppeigensolver.addVectorField(observables[name])
    def run(self):
        self.cppeigensolver.getAllEigenPairs(self.package)
    
    def update(self,frame_num):
        self.pyhessian.update(frame_num)
        self.cppeigensolver.setHessian(self.pyhessian.cpphessian)

class hessian(object):
    def __init__(self,package):
        #Initialize system data and pair potential of the system
        self.package = package;
        self.cppmanager = _elasticitypetsc.PETScManager();
        dimensions = pyglasstools.get_sysdata().pysimbox.dim
        if (package == "petsc"):
            if dimensions == 2:
                self.cpphessian = _elasticitypetsc.PETScHessian2D(pyglasstools.get_sysdata().cppparticledata,pyglasstools.get_potential().cpppairpotential,self.cppmanager,comm)
            else:
                self.cpphessian = _elasticitypetsc.PETScHessian3D(pyglasstools.get_sysdata().cppparticledata,pyglasstools.get_potential().cpppairpotential,self.cppmanager,comm)
        elif (package == "slepc"):
            if dimensions == 2:
                self.cpphessian = _elasticitypetsc.SLEPcHessian2D(pyglasstools.get_sysdata().cppparticledata,pyglasstools.get_potential().cpppairpotential,self.cppmanager,comm)
            else:
                self.cpphessian = _elasticitypetsc.SLEPcHessian3D(pyglasstools.get_sysdata().cppparticledata,pyglasstools.get_potential().cpppairpotential,self.cppmanager,comm)
        self.frame_num = pyglasstools.get_sysdata().frame_num
   
    def update(self,frame_num):
        pyglasstools.get_sysdata().update(frame_num);
        self.cpphessian.setSystemData(pyglasstools.get_sysdata().cppparticledata)
        self.cpphessian.destroyObjects()
        self.cpphessian.assembleObjects()
        self.frame_num = frame_num

        #Another way to update the hessian, by deleteing the object and creating it again from scratch
        #del self.cpphessian
        #if (self.package == "petsc"):
        #    self.cpphessian = _elasticitypetsc.PETScHessian2D(pyglasstools.get_sysdata().cppparticledata,self.pyglasstools.get_potential().cpppairpotential,self.cppmanager,comm)
        #elif (self.package == "slepc"):
        #    self.cpphessian = _elasticitypetsc.SLEPcHessian2D(pyglasstools.get_sysdata().cppparticledata,self.pyglasstools.get_potential().cpppairpotential,self.cppmanager,comm)
        #if self.frame_num != frame_num:
        #    self.frame_num = frame_num
    
    def check_diagonals(self):
        return self.cpphessian.areDiagonalsNonZero()

########
#For I/O Related Routines

class logfile(object):
    global loggers_list

    def __init__(self, filename, names, solver):
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
        self.global_obs = initialize_global(self.__obs_names, pyglasstools.get_sysdata().pysimbox.dim)
        self.solver.add_observables(self.global_obs)
        
        Dim = pyglasstools.get_sysdata().pysimbox.dim
        
        #Create column headers
        if rank == 0 and (not os.path.exists(pyglasstools.get_sysdata().checkpointfile) or (os.path.exists(pyglasstools.get_sysdata().checkpointfile) and os.path.getsize(pyglasstools.get_sysdata().checkpointfile) == 0)):

            self.file.write_shared("{} ".format("Frame"))
            if not (not self.__obs_names):
                for name in self.__obs_names:
                    self.file.write_shared("{} ".format(name))
            self.file.write_shared("\n")
         
    def save(self,frame_num):
        if rank == 0 and self.solver.pyhessian.check_diagonals():
            self.file.write_shared("{} ".format(frame_num))
            for name in self.__obs_names:
                Dim = pyglasstools.get_sysdata().pysimbox.dim
                if "eigenvalue" in name or "nconv" in name or "eigenrelerr" in name:
                    self.global_obs[name].save(self.file)
                if "nonaffinetensor" in name:
                    i = int(name[-4]) 
                    j = int(name[-3])
                    k = int(name[-2]) 
                    l = int(name[-1])
                    self.global_obs['nonaffinetensor'].save(self.file,l+Dim*(k+Dim*(j+Dim*i)))
                self.file.write_shared(" ")
            self.file.write_shared("\n")
    
class fieldlogger(object):
    global loggers_list

    def __init__(self, keyword, names, solver):
        #Save filename
        self.keyword = keyword
        #Next, parse the list of names based on what type of obsercables they are
        self.__obs_names = names
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
        self.field_obs = initialize_field(self.__obs_names, pyglasstools.get_sysdata().pysimbox.dim)
        self.solver.add_observables(self.field_obs)
        self.file = _pyglasstools.MPILogFile(comm, "{}".format(keyword)+".xyz")
        
        #Now, we specify the size of particles assigned to a process
        #ends = [sum(counts[:p+1]) for p in range(nprocs)]

    def save_perrank(self,frame_num):
        Dim = pyglasstools.get_sysdata().pysimbox.dim

        ave, res = divmod(len(pyglasstools.get_sysdata().traj[frame_num].particles.position), size)
        counts = [ave + 1 if p < res else ave for p in range(size)]
        # determine the starting and ending indices of each sub-task
        starts = [sum(counts[:p]) for p in range(size)]
        ends = [sum(counts[:p+1]) for p in range(size)]

        # converts data into a list of arrays 
        #Do scatter for each held eigenvectors
        for name in self.__obs_names:
            if "eigenvector" in name or "displacement" in name or "loclandscape" in name:
                self.field_obs[name].scatterVector(starts[rank],counts[rank])
        for i in range(starts[rank],ends[rank]):
            outline = "{} ".format(1)
            outline += "{} ".format(pyglasstools.get_sysdata().traj[frame_num].particles.position[i,0])
            outline += "{} ".format(pyglasstools.get_sysdata().traj[frame_num].particles.position[i,1])
            if Dim == 3:
                outline += "{} ".format(pyglasstools.get_sysdata().traj[frame_num].particles.position[i,2])
            outline += "{} ".format(0.5*pyglasstools.get_sysdata().traj[frame_num].particles.diameter[i])
            for name in self.__obs_names:
                flattenindex = 0;
                if "displacement" in name or "loclandscape" in name:
                    outline += self.field_obs[name].gettostring(i-starts[rank])
                if "eigenvector" in name:
                    outline += self.field_obs[name].gettostring(i-starts[rank])
            outline += "\n"
            self.file.write_shared(outline);
    
    def get_vectorfield(self, name,frame_num):
        vector = []
        ave, res = divmod(len(pyglasstools.get_sysdata().traj[frame_num].particles.position), size)
        counts = [ave + 1 if p < res else ave for p in range(size)]
        starts = sum(counts[:rank])
        ends = sum(counts[:rank+1])
        for i in range(starts,ends):
            vector.append(np.reshape(self.field_obs[name].getVectorValue(i),(3,)))
        #Next, we gather all of these vectors
        newvector = comm.all_gather_v(vector)
        vector = [item for sublist in newvector for item in sublist]
        #vector = np.reshape(newvector,(len(pyglasstools.get_sysdata().traj[frame_num].particles.position),3));
        del newvector;
        return vector#Reshape
        
    def save(self,frame_num):
        if self.solver.pyhessian.check_diagonals():
            if rank == 0:
                self.file.write_shared("{:d} \n".format(len(pyglasstools.get_sysdata().traj[frame_num].particles.position)))
                self.file.write_shared("#Frame {:d}  \n".format(frame_num))
            comm.barrier()
            self.save_perrank(frame_num)
