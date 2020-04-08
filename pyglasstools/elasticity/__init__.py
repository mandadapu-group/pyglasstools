from pyglasstools.elasticity import _elasticity
import numpy as np

### Try to import petsc4py and slepc4py
isslepc = True;
try:
    import sys, slepc4py
    slepc4py.init(sys.argv)
except ImportError:
    print("[WARNING] No slepc4py installation found. eigs_slepc will be unusable")
    isslepc = False;

ispetsc = True;
try:
    from petsc4py import PETSc
except ImportError:
    print("[WARNING] No petsc4py installation found. eigs_slepc will be unusable")
    ispetsc = False;


class hessian(object):
    def __init__(self, sysdata,potential):
        self.H = _elasticity.Hessian(sysdata._getParticleSystem(),potential._getPairPotential())
    
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
    
    def eigs(self, selrule = 'LM', nev = 1, ncv = 5, maxiter=1000, tol=1e-10):
        self.H.getEigenDecomposition(selrule,nev,ncv,maxiter,tol);
    
    def eigs_slepc(self, maxiter=1000, tol=1e-10):
        if (isslepc == True and ispetsc ==True):
            SLEPc = slepc4py.SLEPc
            Print = PETSc.Sys.Print
            
            #Form PETSc Matrix
            Print("Build PETSc Sparse matrix")
            A = PETSc.Mat()
            A.createAIJ(size=self.H.hessian.shape, csr=(self.H.hessian.indptr, self.H.hessian.indices, self.H.hessian.data),comm=PETSc.COMM_WORLD) 
            A.setUp(); A.assemble()
            A.setType('seqaij')
            opts = PETSc.Options()
            
            Print("Build SLEPc Eigensolver")
            E = SLEPc.EPS(); E.create()
            E.setOperators(A)
            E.setProblemType(SLEPc.EPS.ProblemType.HEP)
            E.setFromOptions()
            E.setTolerances(tol,maxiter) 
            
            xdir, ydir = A.getVecs()
            
            if any("eps_interval" in s for s in sys.argv):
                Print("Peforming spectrum slicing. No deflation of null space is needed")
            else: 
                Print("Deflate Null space")
                for i in range(self.H.hessian.shape[0]):
                    if (i % 2):
                        xdir[i] = 1
                        ydir[i] = 0
                    else:
                        xdir[i] = 0
                        ydir[i] = 1
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
            Print("Number of converged eigenpairs %d" % nconv)
            if nconv > 0:
                # Create the results vectors
                vr, wr = A.getVecs()
                vi, wi = A.getVecs()
                #
                Print()
                Print("        k          ||Ax-kx||/||kx|| ")
                Print("----------------- ------------------")

                tempeigvecs = []
                tempeigvals = []
                for i in range(nconv):
                    k = E.getEigenpair(i, vr, vi)
                    tempeigvecs.append(vr[:])
                    tempeigvals.append(k)
                self.H.eigenvecs = np.array(tempeigvecs).T
                self.H.eigenvals = np.array(tempeigvals)
        else:
            print("[WARNING] petsc4py and slepc4py are not installed. This function won't do anything")

    def _getObservable(self):
        return self.H
