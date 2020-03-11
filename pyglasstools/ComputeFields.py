import numpy as np
import numba as nb

@nb.njit(nogil=True,cache=True)
def PBC(dx,box):
    for i in range(dx.shape[0]):
        if (dx[i] > box*0.5):
            dx[i] -= box
        elif (dx[i] <= -box*0.5):
            dx[i] += box
    return dx

from scipy.spatial import cKDTree

@nb.njit(nogil=True,fastmath=True,cache=True)
def Delta(x,rprime,rcut,boxsize):
    dr = PBC(x-rprime,boxsize)#/rcut
    absr = np.sqrt(np.sum(dr**2))/rcut#-rcut
    if absr < 1:
        return 15/(8*np.pi*rcut**2)*(1-2*absr**4+absr**8)
    else:
        return 0

@nb.njit(nogil=True,fastmath=True,cache=True)
def bond(s,x,ri,dr,rcut,boxsize):
    return Delta(x,ri+s*dr,rcut,boxsize)

weightquad = np.array([0.568888889,0.47862867,0.47862867,0.236926885,0.236926885])
xquad = np.array([0,-0.53846931,0.53846931,-0.906179846, 0.906179846])     
@nb.njit(nogil=True,fastmath=True,cache=True)
def integrate_bond(x,ri,dr,rcut,boxsize):
    y = (xquad+1)/2.0 
    val = 0
    for i in range(5):
        val += 0.5*(weightquad[i]*bond(y[i],x,ri,dr,rcut,boxsize))
    return val


def getPairArrays(rtree, grids,boxsize,radius):
    pair = rtree.query_ball_point(grids,radius)
    pairinf = []
    newpair = []
    count = 0
    for i in range(len(pair)):
        pairinf.append([i+count,len(pair[i])])
        for item in pair[i]:
            newpair.append(item)
        count += len(pair[i])-1;
        #print(count)
    pair = np.array(newpair).flatten('C')
    pairinf = np.array(pairinf)
    return pair, pairinf

def getCGField(cg_rcut,force_rcut,r,diameters,boxsize,dx,compute_force,forceargs=()):
    #Construct Tree
    rtree = cKDTree(r, boxsize=boxsize)
    #Find Pairs Pertaining to Force Computations 
    forcepair, forcepair_inf = getPairArrays(rtree, r,boxsize,force_rcut)    
    #Compute Grid Points and Find Pairs Pertaining to Coarse Graining Function
    nmax = int(boxsize/dx)
    points = np.linspace(-boxsize/2.0,+boxsize/2.0,nmax)#/boxsize
    x = []
    for i in range(len(points)): 
        for j in range(len(points)):
            x.append([points[i],points[j]])
    x = np.asarray(x,dtype=np.float32)
    
    #Find Pairs Pertaining to Force Computations 
    cgpair, cgpair_inf = getPairArrays(rtree, x,boxsize,cg_rcut)    
    
    xsize = len(x)
    rho = np.zeros(xsize,dtype=float)
    Tv = np.zeros((xsize,3),dtype=float)
    rho, Tv = computeCGField(rho,Tv,r,diameters,forcepair,forcepair_inf,cgpair,cgpair_inf,x,xsize,boxsize,force_rcut,cg_rcut,compute_force,*forceargs)
    return rho, Tv    

#@nb.njit(nogil=True,fastmath=True,parallel=True)
@nb.njit(nogil=True,fastmath=True)
def computeCGField(rho,Tv,r,diameters,forcepair,forcepair_inf,cgpair,cgpair_inf,x,xsize,boxsize,force_rcut,cg_rcut,compute_force,*forceargs):
    #for n in nb.prange(xsize):
    for n in range(xsize):
        pair_ni = cgpair[cgpair_inf[n,0]:np.sum(cgpair_inf[n])]
        for i in pair_ni:
            rho[n] += Delta(x[n],r[i],cg_rcut,boxsize)
            pair_ij = forcepair[forcepair_inf[i,0]:np.sum(forcepair_inf[i])]
            for j in pair_ij:
                dr = PBC(r[j]-r[i],boxsize) 
                if (dr[0] != 0 and dr[1] != 0):
                    F = compute_force(dr,diameters[i],diameters[j],*forceargs) 
                    if F[0] != 0 or F[1] != 0:
                        prefactor = integrate_bond(x[n],r[i],dr,cg_rcut,boxsize)
                        if prefactor == 0 or prefactor < 1e-6:
                            continue
                        else:
                            Tv[n,0] += 0.5*F[0]*dr[1]*prefactor
                            Tv[n,1] += 0.5*dr[0]*F[0]*prefactor
                            Tv[n,2] += 0.5*dr[1]*F[1]*prefactor
                    else:
                        continue
                else:
                    continue
    return rho, Tv
