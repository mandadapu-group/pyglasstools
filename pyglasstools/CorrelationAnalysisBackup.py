import numpy as np
import numba as nb
import time
from tqdm import tqdm
import os
import pyglasstools.BlockAverage as block
#import BlockAverage as block

@nb.njit(nogil=True,cache=True)
def PBC(dx,box):
    for i in range(dx.shape[0]):
        if (dx[i] > box*0.5):
            dx[i] -= box
        elif (dx[i] <= -box*0.5):
            dx[i] += box
    return dx

@nb.njit(nogil=True, fastmath=True)
def getCorrelAtTimeT(positions,datapoints,index,ordparam,args=()):
    val = np.zeros(datapoints-index)
    for i in nb.prange(datapoints-index):
        pos0 = positions[i]
        post = positions[i+index]
        val[i] = ordparam(post,pos0,*args)
    return np.mean(val);#np.mean(val)

#@nb.njit(nogil=True, fastmath=True)
def getCorrelFunc(positions,datapoints,ordparam,args=()):
    x = ordparam(positions,datapoints,*args)
    xp = ifftshift((x - np.average(x))/np.std(x))
    n, = xp.shape
    xp = np.r_[xp[:n//2], np.zeros_like(xp), xp[n//2:]]
    f = fft(xp)
    p = np.absolute(f)**2
    pi = ifft(p)
    return np.real(pi)[:n//2]/(np.arange(n//2)[::-1]+n//2)
    #for i in nb.prange(datapoints-index):
    #    pos0 = positions[i]
    #    post = positions[i+index]
    #    val[i] = ordparam(post,pos0,*args)
    #return np.mean(val);#np.mean(val)

#def getCorrelAtTimeT(positions,datapoints,index,ordparam,args=()):
#    val = getCorrelSample(positions,datapoints,index,ordparam,*args)
    
    #Perform Block Error Analysis
#    blockerr = block.BlockAverage(val)
#    (n, num, err, err_err) = blockerr.get_hierarchical_errors()
#    (index, err_est) = blockerr.get_error_estimate()
    #Return Mean and Error Bar
#    return np.mean(val)#, err_est 
@nb.njit(nogil=True, fastmath=True,parallel=True)
def densalt(pos,datapoints,wavek,k_sample,n_particles,boxsize):
    rho = np.zeros(datapoints)#,dtype=complex)
    for k in range(datapoints):
        print(k)
        realpart = np.zeros(n_particles)#,dtype=complex)
        for i in nb.prange(n_particles):
            for j in range(k_sample):
                dr = pos[k,i];#PBC(post[i]-pos0[i],boxsize)
                realpart[i] +=  np.cos(np.dot(wavek[j],dr))#+1j*np.sin(np.dot(wavek[j],dr))
        rho[k] = np.mean(realpart)/k_sample
    return rho
@nb.njit(nogil=True, fastmath=True,parallel=True)
def density(post,pos0,k,k_sample,n_particles,boxsize):
    realpart = np.zeros(n_particles)
    for i in nb.prange(n_particles):
        for j in range(k_sample):
            dr = PBC(post[i]-pos0[i],boxsize)
            realpart[i] +=  np.cos(np.dot(k[j],dr))
    return np.mean(realpart)/k_sample

def getISF(trajpos,snapshots,kpeak,k_sample,n_particles,boxsize,mode="write"):
    theta = np.linspace(0,2*np.pi-1e-3,k_sample)
    wavek = np.zeros((k_sample,3),dtype=np.float32)
    for i in range(k_sample): 
        wavek[i,0] = kpeak*np.cos(theta[i])
        wavek[i,1] = kpeak*np.sin(theta[i])
        wavek[i,2] = 0
    if (mode == "write"):
        with open("fsk.txt","w") as file:
            fsk = getCorrelFunc(trajpos,snapshots,densalt,args=(wavek,k_sample,n_particles,boxsize))
            fsk[0] = 1.0
            file.write("Fs(k,t) with k = {} \n".format(kpeak))
            file.write("Frame Val \n")
            file.write("{} {} \n".format(0, fsk[0])) 
            for i in tqdm(range(1,snapshots), desc="Computing Fs(k,t)"):
                #fsk[i] = getCorrelAtTimeT(trajpos,snapshots,i,density,args=(wavek,k_sample,n_particles,boxsize))
                file.write("{} {} \n".format(i, fsk[i]))
                file.flush()
                os.fsync(file.fileno())
        return fsk
    else:
        fsk[0] = 1.0
        #fsk_err[0] = 0
        for i in tqdm(range(1,snapshots), desc="Computing Fs(k,t)"):
            fsk[i] = getCorrelAtTimeT(trajpos,snapshots,i,density,args=(wavek,k_sample,n_particles,boxsize))
        return fsk
