import numpy as np
import numba as nb
import time
from tqdm import tqdm
import os
import pyglasstools.BlockAverage as block
import pyglasstools.ComputeFields as field
from scipy.spatial import cKDTree
#import BlockAverage as block

import pyfftw
def autocorrelation_fft(x):
    a = np.concatenate((x,np.zeros(x.shape[0]-1))) 
    A = pyfftw.interfaces.numpy_fft.fft(a,threads=1)
    S = A*np.conj(A)
    c_fourier = pyfftw.interfaces.numpy_fft.ifft(S,threads=1)
    c_fourier = c_fourier[:(c_fourier.size//2)+1]
    return np.real(c_fourier)/(np.arange(x.size, 0, -1) )

#@nb.jit(nogil=True, fastmath=True)#,forceobj=True)
@nb.njit(nogil=True,fastmath=True)#,forceobj=True)
def autocorrelation_timeavg(x,xconj,N):
    val = np.zeros_like(xconj)
    for r in range(N):
        for k in range(N-r):
            val[r] += x[k]*xconj[k+r]
        val[r] /= float(N-r)
    return val

#from scipy.special import airy
@nb.jit(nogil=True, fastmath=True,forceobj=True)
def getLocalCorrelFunc(positions,datapoints,n_particles,wavek,k_sample,boxsize,ordparam,mode,args=()):
    val = np.zeros((n_particles,datapoints),dtype=complex)
    for i in tqdm(range(n_particles), desc="Computing Auto-Correlation Function : "):
        for j in range(datapoints):
            val[i,j] = np.sum(ordparam(i,positions[j,:,:],wavek,*args))
        if (mode == "fft"):
            val[i,:] = autocorrelation_fft(val[i,:])
        elif (mode == "timeavg"):
            valconj = np.conj(val[i,:])
            val[i,:] = autocorrelation_timeavg(val[i,:],valconj,datapoints)
        else:
            raise ValueError('mode for correlation function calculation is undefined. Only ''fft'' and ''timeavg'' are available.')
    return np.mean(val/k_sample,axis=0);#np.mean(val)
"""
#@nb.jit(nogil=True, fastmath=True,forceobj=True)
def getGlobalPairCorrelFunc(positions,datapoints,n_particles,wavek,k_sample,ordparam,mode,args=()):
    val = np.zeros(datapoints,dtype=complex)
    for i in tqdm(range(datapoints), desc="Computing Global Observable : "):
        for j in range(k_sample):
            val[i] += np.sum(ordparam(positions[i,:,:].T,wavek[j],*args))
    print(val)
    #val /= k_sample
    if (mode == "fft"):
        val = autocorrelation_fft(val)
    elif (mode == "timeavg"):
        valconj = np.conj(val)
        val = autocorrelation_timeavg(val,valconj,datapoints)
    else:
        raise ValueError('mode for correlation function calculation is undefined. Only ''fft'' and ''timeavg'' are available.')
    return val/(k_sample*n_particles);#np.mean(val)
#from scipy.special import airy
@nb.jit(nogil=True, fastmath=True,forceobj=True)
def getGlobalPairCorrelFunc(positions,datapoints,n_particles,wavek,k_sample,trajpair,trajpair_inf,ordparam,mode,args=()):
    val = np.zeros((n_particles,datapoints),dtype=complex)
    for i in tqdm(range(datapoints), desc="Computing Global Observable : "):
        for j in range(k_sample):
            val[i,j] = np.sum(ordparam(i,positions[i,:,:],wavek,x,trajpair[i],trajpair_inf[i],*args))
    #for i in tqdm(range(n_particles), desc="Computing Auto-Correlation Function : "):
        #x = np.zeros(k_sample,dtype=complex)
        #for j in range(datapoints):
            #Construct Tree
        if (mode == "fft"):
            val[i,:] = autocorrelation_fft(val[i,:])
        elif (mode == "timeavg"):
            valconj = np.conj(val[i,:])
            val[i,:] = autocorrelation_timeavg(val[i,:],valconj,datapoints)
        else:
            raise ValueError('mode for correlation function calculation is undefined. Only ''fft'' and ''timeavg'' are available.')
    return np.mean(val/k_sample,axis=0);#np.mean(val)
"""
@nb.njit(nogil=True, fastmath=True)
def correctPositions(positions,n_particles,datapoints,boxsize):
    for i in range(n_particles):
        for j in range(1,datapoints):
            dx = positions[j,i,:]-positions[j-1,i,:]
            for k in range(dx.shape[0]):
                if (dx[k] > boxsize*0.5):
                    positions[j,i,k] -= boxsize
                elif (dx[k] <= -boxsize*0.5):
                    positions[j-1,i,k] += boxsize
    return positions

@nb.njit(nogil=True, fastmath=True,parallel=True)
def localdensity(index,pos,wavek):
    return np.exp(1j*np.dot(wavek,pos[index]))
"""
@nb.njit(nogil=True, fastmath=True)
def localstress(index,pos,wavek,x,forcepair,forcepair_inf,boxsize,component,diameters,compute_force,*forceargs):
    pair_ij = forcepair[forcepair_inf[index,0]:np.sum(forcepair_inf[index])]
    for j in pair_ij:
        dr = field.PBC(pos[j]-pos[index],boxsize) 
        if (dr[0] != 0 and dr[1] != 0):
            F = compute_force(dr,diameters[index],diameters[j],*forceargs) 
            if F[0] > 0 and F[1] > 0:
                kdottedr = np.dot(wavek,dr)
                a = 0.5*np.fmod(dr,boxsize*0.5)+pos[index]
                kdottedmid = np.dot(wavek,a)
                prefactor = np.exp(-1j*kdottedmid)*np.sin(0.5*kdottedr)/(kdottedr)#np.exp(#field.integrate_bond(pos[index],pos[index],dr,cg_rcut,boxsize)
                if (component == 0):
                    x += F[0]*dr[0]*prefactor
                elif (component == 1):
                    x +=  F[1]*dr[1]*prefactor
                elif (component == 2):
                    x += F[0]*dr[1]*prefactor
            else:
                continue
        else:
            continue
    return x
def getSelfStressF(traj,kpeak,k_sample,component,compute_force,forceargs=(),mode="fft",filename="selfstressAF.txt"):
    #Turn it into proper array
    snapshots = len(traj)
    n_particles = traj[0].particles.N
    boxsize = np.float64(traj[0].configuration.box[0])
    diameters = traj[0].particles.diameter
    trajpos = []
    trajforcepair = []
    trajforcepair_inf = []
    force_rcut = 2**(1/6.0)*1.4
    for i in tqdm(range(snapshots),desc="Generating Neighbor List for Trajectory"):
        r = traj[i].particles.position
        trajpos.append(r.astype(np.float64))
        rtree = cKDTree(r[:,0:2]+0.5*boxsize, boxsize=boxsize)
        #Find Pairs Pertaining to Force Computations 
        forcepair, forcepair_inf = field.getPairArrays(rtree,r[:,0:2]+0.5*boxsize,boxsize,force_rcut)
        trajforcepair.append(forcepair)
        trajforcepair_inf.append(forcepair_inf)
    trajpos = np.asarray(trajpos,dtype=np.float64)
    

    theta = np.linspace(0,2*np.pi-1e-3,k_sample)
    wavek = np.zeros((k_sample,3),dtype=np.float64)
    for i in range(k_sample): 
        wavek[i,0] = kpeak*np.cos(theta[i])
        wavek[i,1] = kpeak*np.sin(theta[i])
        wavek[i,2] = 0
    
    #print("Computing Particle Positions corrected by Periodic Boundary Conditions")
    #trajpos =  correctPositions(trajpos,n_particles,snapshots,boxsize)
    #rtree = cKDTree(r, boxsize=boxsize)
    stressaf = getLocalPairCorrelFunc(trajpos,snapshots,n_particles,wavek,k_sample,trajforcepair,trajforcepair_inf,localstress,mode,args=(boxsize,component,diameters,compute_force,*forceargs))
    stressaf /=  stressaf[0]
    print("Writing Output to {}".format(filename))
    with open("{}".format(filename),"w") as file:
        file.write("Fs(k,t) with k = {} \n".format(kpeak))
        file.write("Frame Val \n")
        for i in range(snapshots):
            file.write("{} {} \n".format(i, np.abs(stressaf[i])))
            file.flush()
            os.fsync(file.fileno())
    return stressaf
"""

def getISF(traj,kpeak,k_sample,mode="fft",filename="fsk.txt"):
    snapshots = len(traj)
    n_particles = traj[0].particles.N
    boxsize = traj[0].configuration.box[0]
    #Turn it into proper array
    trajpos = []
    for i in range(len(traj)):
        r = traj[i].particles.position#+0.5*boxsize
        trajpos.append(r)
    
    trajpos = np.asarray(trajpos)
    
    
    theta = np.linspace(0,2*np.pi-1e-3,k_sample)
    wavek = np.zeros((k_sample,3),dtype=np.float32)
    for i in range(k_sample): 
        wavek[i,0] = kpeak*np.cos(theta[i])
        wavek[i,1] = kpeak*np.sin(theta[i])
        wavek[i,2] = 0
    
    fsk = getLocalCorrelFunc(trajpos,snapshots,n_particles,wavek,k_sample,boxsize,localdensity,mode)
    
    print("Writing Output to {}".format(filename))
    with open("{}".format(filename),"w") as file:
        file.write("Fs(k,t) with k = {} \n".format(kpeak))
        file.write("Frame Val \n")
        for i in range(snapshots):
            file.write("{} {} \n".format(i, np.abs(fsk[i])))
            file.flush()
            os.fsync(file.fileno())
    return fsk


def checkAgingISF(trajpos,snapshots,kpeak,k_sample,n_particles,boxsize,numoffsk,mode="fft",filename="fsk_aging.txt"):
    interval = int(snapshots//numoffsk)
    fsk = []
    count = 0
    for j in range(numoffsk):
        print("Initial Frame: {}".format(j*interval))
        if(j*interval > snapshots // 2):
            print("Number of Frames is Less than Half of Your Original Trajectory. Exiting Aging Calculation")
            break
        else:
            fsk.append(getISF(trajpos[j*interval:],snapshots-j*interval,kpeak,k_sample,n_particles,boxsize,mode,filename="{}.{}".format(filename,j)))
    return fsk
