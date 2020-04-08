from pyglasstools.elasticity import _elasticity

class hessian(object):
    def __init__(self, sysdata,potential):
        self.H = _elasticity.Hessian(sysdata._getParticleSystem(),potential._getPairPotential())
    #Redefine attributes so that it directly access SimBox C++ class 
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
        return self.H.getEigenDecomposition(selrule,nev,ncv,maxiter,tol);

    def _getObservable(self):
        return self.H
