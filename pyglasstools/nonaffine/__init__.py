from pyglasstools.nonaffine import _nonaffine
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

### Try to import petsc4py and slepc4py
isslepc = True;
try:
    import sys, slepc4py
    slepc4py.init(sys.argv)
except ImportError:
    if rank == 0:
        print("[WARNING] No slepc4py installation found. eigs_slepc will be unusable")
    isslepc = False;

ispetsc = True;
try:
    from petsc4py import PETSc
    Print = PETSc.Sys.Print
except ImportError:
    if rank == 0:
        print("[WARNING] No petsc4py installation found. eigs_slepc will be unusable")
    ispetsc = False;

class hessian(object):
    def __init__(self, sysdata,potential):
        if rank > 0:
            print("[WARNING] Performing MPI run. The hessian class supports no MPI parallelization.")
            print("[WARNING] MPI parallelization only exists in hessian_slepc class.")
            if (ispetsc == False or isslepc == False):
                print("[WARNING] Is petsc4py installed? {} Is slepc4py installed? {}".format(ispetsc,isslepc))
            elif (ispetsc == True and isslepc == True):
                print("[WARNING] petsc4py and slepc4py are installed. Consider using hessian_slepc class instead".format(ispetsc,isslepc))
        self.H = _nonaffine.Hessian(sysdata._getParticleSystem(),potential._getPairPotential())
     
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
        if rank == 0:
            self.H.getEigenDecomposition(selrule,nev,ncv,maxiter,tol);
    
    ## Shortcut Function to Compute All Eigenpairs
    def alleigs(self,dim,maxiter=10000, tol=1e-8):
        if rank == 0:
            self.H.getEigenDecomposition('LM',self.H.hessian.shape[0]-dim,self.H.hessian.shape[0],maxiter,tol);
    
    ## Building pseudoinverse
    def build_pinv(self,tol):
        if rank == 0:
            #Check the number of eigenpairs, and compute the Frobenius norm error
            self.H.buildPseudoInverse(tol)
    
    ## Check if system full eigendecomposition reproduces the Hessian
    def check_alleigs(self):
        if rank == 0:
            self.H.checkFullDecompError()

    def _getObservable(self):
        return self.H


##Helper function to convert scipy sparse matrix into a PETSc matrix
def csrmatrix2PETScMat(L):
    """
    Converts a sequential scipy sparse matrix (on process 0) to a PETSc
    Mat ('aij') matrix distributed on all processes
    input : L, scipy sparse matrix on proc 0
    output: PETSc matrix distributed on all procs
    """

    # Get the data from the sequential scipy matrix
    if rank == 0:
        if L.format == 'csr':
            L2 = L
        else:
            L2 = L.tocsr()
        Ai  = L2.indptr
        Aj  = L2.indices
        Av  = L2.data
        nnz = len(Aj)
        n,m = L2.shape
    else:
        n   = None
        m   = None
        nnz = None
        Ai  = None
        Aj  = None
        Av  = None

    # Broadcast sizes
    n   = comm.bcast(n  ,root = 0)
    m   = comm.bcast(m  ,root = 0)
    nnz = comm.bcast(nnz,root = 0)

    B = PETSc.Mat()
    B.create(comm)
    B.setSizes([n, m])
    B.setType('aij')
    B.setFromOptions()

    # Create a vector to get the local sizes, so that preallocation can be done later
    V = PETSc.Vec()
    V.create(comm)
    V.setSizes(n)
    V.setFromOptions()
    istart,iend = V.getOwnershipRange()
    V.destroy()

    nloc = iend - istart

    Istart = comm.gather(istart,root = 0)
    Iend   = comm.gather(iend  ,root = 0)

    if rank == 0:
        nnzloc = np.zeros(comm.size,'int')
        for i in range(comm.size):
            j0        = Ai[Istart[i]]
            j1        = Ai[Iend  [i]]
            nnzloc[i] = j1 - j0
    else:
        nnzloc = None

    nnzloc = comm.scatter(nnzloc,root = 0)

    ai = np.zeros(nloc+1   ,PETSc.IntType)
    aj = np.zeros(nnzloc+1 ,PETSc.IntType)
    av = np.zeros(nnzloc+1 ,PETSc.ScalarType)

    if rank == 0:
        j0        = Ai[Istart[0]]
        j1        = Ai[Iend  [0]]
        ai[:nloc  ] = Ai[:nloc]
        aj[:nnzloc] = Aj[j0:j1]
        av[:nnzloc] = Av[j0:j1]

    for iproc in range(1,comm.size):
        if rank == 0:
            i0        = Istart[iproc]
            i1        = Iend  [iproc]
            j0        = Ai[i0]
            j1        = Ai[i1]
            comm.Send(Ai[i0:i1], dest=iproc, tag=77)
            comm.Send(Aj[j0:j1], dest=iproc, tag=78)
            comm.Send(Av[j0:j1], dest=iproc, tag=79)
        elif rank == iproc:
            comm.Recv(ai[:nloc  ], source=0, tag=77)
            comm.Recv(aj[:nnzloc], source=0, tag=78)
            comm.Recv(av[:nnzloc], source=0, tag=79)

    ai = ai- ai[0]
    ai[-1] = nnzloc+1

    B.setPreallocationCSR((ai,aj))
    B.setValuesCSR(ai,aj,av)
    B.assemble()

    return B


class hessian_slepc(object):
    def __init__(self, sysdata,potential):
        self.H = _nonaffine.Hessian(sysdata._getParticleSystem(),potential._getPairPotential())
    
    #Redefine attributes so that it directly access Hessian C++ class 
    #attributes
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
    
    ## Building pseudoinverse
    def build_pinv(self,tol):
        if rank == 0:
            #Check the number of eigenpairs, and compute the Frobenius norm error
            self.H.buildPseudoInverse(tol)
    ## Implementation of eigendecomposition using SLEPc and PETSc
    def eigs(self, maxiter=1000, tol=1e-10):
        if (isslepc == True and ispetsc ==True):
            SLEPc = slepc4py.SLEPc
            
            #Form PETSc Matrix
            Print("Build PETSc Sparse matrix")
            A = csrmatrix2PETScMat(self.H.hessian)#PETSc.Mat()
            
            Print("Build SLEPc Eigensolver")
            E = SLEPc.EPS(); E.create()
            E.setOperators(A)
            E.setProblemType(SLEPc.EPS.ProblemType.HEP)
            E.setFromOptions()
            E.setTolerances(tol,maxiter) 
            
            xdir, ydir = A.getVecs()
            istart, iend = xdir.getOwnershipRange() 
            if any("eps_interval" in s for s in sys.argv):
                Print("Peforming spectrum slicing. No deflation of null space is needed")
            else: 
                Print("Deflate Null Space")
                for i in range(istart,iend):
                    if (i % 2):
                        xdir[i] = 1
                        ydir[i] = 0
                    else:
                        xdir[i] = 0
                        ydir[i] = 1
                xdir.assemble()
                ydir.assemble()
                E.setDeflationSpace([xdir,ydir])
            
            Print("Solve!")
            E.solve()
            Print()
            Print("******************************")
            Print("*** SLEPc Solution Results ***")
            Print("******************************")
            Print()

            its = E.getIterationNumber()
            Print("Number of iterations of the method: %d" % its)

            eps_type = E.getType()
            Print("Solution method: %s" % eps_type)

            nev, ncv, mpd = E.getDimensions()
            Print("Number of requested eigenvalues: %d" % nev)

            tol, maxit = E.getTolerances()
            Print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))

            nconv = E.getConverged()
            self.H.nconv = nconv
            Print("Number of converged eigenpairs %d" % nconv)
            if nconv > 0:
                # Create the results vectors
                vr, wr = A.getVecs()
                vi, wi = A.getVecs()
                #
                #Print()
                #Print("        k          ||Ax-kx||/||kx|| ")
                #Print("----------------- ------------------")
                
                istart, iend = vr.getOwnershipRange()
                if rank == 0:
                    tempeigvecs = []
                    tempeigvals = []
                else:
                    tempeigvecs = None
                    tempeigvals = None

                for i in range(nconv):
                    k = E.getEigenpair(i, vr, vi)
                    error = E.computeError(i)
                  #  if k.imag != 0.0:
                  #      Print(" %9f%+9f j %12g" % (k.real, k.imag, error))
                  #  else:
                  #      Print(" %12f      %12g" % (k.real, error))
                    Vr = comm.gather(vr[istart:iend],root = 0)
                    if rank == 0:
                        Vr = np.concatenate(Vr)
                        tempeigvecs.append(Vr[:])
                        tempeigvals.append(np.real(k))
                
                #Destroy Eigensolver and Matrix object because data is already stored onto Root Process
                E.destroy()
                A.destroy()

                #I can then handle I/O here. Perhaps there is a better way . . . .                 
                Print("Copying Eigenvectors and Eigenvalues to Root Process")
                if rank == 0:
                    self.H.eigenvecs = np.array(tempeigvecs).T
                    self.H.eigenvals = np.array(tempeigvals)
                    del tempeigvecs[:] 
                    del tempeigvals[:]
        else:
            Print("[WARNING] petsc4py and slepc4py are not installed. eigs_slepc will do nothing.")
    
    def check_alleigs(self):
        if rank == 0:
            self.H.checkFullDecompError()
    def _getObservable(self):
        return self.H
