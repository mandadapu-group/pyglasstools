import pyglasstools
from pyglasstools import _pyglasstools, comm, rank, size
import numpy as np
import os
import gsd.hoomd

class simbox(object):
    R""" Define box ndim.

    """
    def __init__(   self, Lx=1.0, Ly=1.0, Lz=1.0,
                    origin=np.zeros(3).astype('float64'), 
                    ndim=3, L=None, volume=None):
        #check dimensions first
        if ndim == 2:
            Lz = 1.0;
        #Initialize the box
        if L is not None:
            self.cppsimbox = _pyglasstools.SimBox(L,origin,ndim);
        else:
            boxsize = np.array([Lx,Ly,Lz]).astype('float64')
            self.cppsimbox = _pyglasstools.SimBox(boxsize,origin,ndim);
    
    #Redefine attributes so that it directly access cppsimbox C++ class 
    #attributes
    def __getattr__(self,attr):
            orig_attr = self.cppsimbox.__getattribute__(attr)
            if callable(orig_attr):
                def hooked(*args, **kwargs):
                    self.pre()
                    result = orig_attr(*args, **kwargs)
                    # prevent cppsimbox from becoming unwrapped
                    if result == self.cppsimbox:
                        return self
                    self.post()
                    return result
                return hooked
            else:
                return orig_attr

def read_gsd(filename,checkpointfile = None,frame_num = 0,ndim=2):
    traj_gsd = gsd.hoomd.open(name=filename,mode='rb')
    return data_gsd(traj_gsd, checkpointfile, frame_num,ndim)

class data_gsd(object):
    global globalsysdata
    R""" Object to store particle data based on GSD format

    """
    def __init__(self, traj_gsd, checkpointfile, frame_num,ndim=2):
        self.frame_num = frame_num
        self.checkpointfile = checkpointfile
        #Store trajectory data
        self.traj = traj_gsd; 
        self.ndim = ndim;
        #Store simulation box
        self.pysimbox = simbox(   Lx=self.traj[frame_num].configuration.box[0],
                                  Ly=self.traj[frame_num].configuration.box[1],
                                  Lz=self.traj[frame_num].configuration.box[2],
                                  ndim=ndim)
        # create the c++ mirror class
        self.cppparticledata = _pyglasstools.ParticleSystem(   self.pysimbox.cppsimbox,
                                                            len(self.traj[frame_num].particles.diameter),
                                                            self.traj[frame_num].particles.diameter,
                                                            self.traj[frame_num].particles.mass,
                                                            self.traj[frame_num].particles.position,
                                                            self.traj[frame_num].particles.velocity 
                                                            );
        pyglasstools.set_sysdata(self)
    
    def setup_checkpoint(self, frame_list):
        #Setup the checkpoint file, if specified 
        if self.checkpointfile is not None and not os.path.exists(self.checkpointfile):
            if rank == 0:
                with open(self.checkpointfile,"w") as f:
                    f.write("Frame {}".format(frame_list[0]))
                    f.close()
        elif self.checkpointfile is not None and os.path.exists(self.checkpointfile):
            with open(self.checkpointfile,"r") as f:
                frame_num = int(f.readline().split(' ')[1])
                frame_list = [i for i in range(frame_num,frame_list[-1]+1)]
                f.close()
        return frame_list

    def update(self,frame_num):
        #Store simulation box
        if frame_num != self.frame_num:
            self.simbox = simbox(   Lx=self.traj[frame_num].configuration.box[0],
                                    Ly=self.traj[frame_num].configuration.box[1],
                                    Lz=self.traj[frame_num].configuration.box[2],
                                    ndim=self.ndim)
            
            # create the c++ mirror class
            self.cppparticledata = _pyglasstools.ParticleSystem(    self.pysimbox.cppsimbox, 
                                                                    len(self.traj[frame_num].particles.diameter),
                                                                    self.traj[frame_num].particles.diameter,
                                                                    self.traj[frame_num].particles.mass,
                                                                    self.traj[frame_num].particles.position,
                                                                    self.traj[frame_num].particles.velocity 
                                                                    );
            self.frame_num = frame_num
    
    def move_particles(self,displacement):
        self.cppparticledata.moveParticles(displacement)
    def move_particles(self):
        self.cppparticledata.moveParticles()
    
    def set_displacement(self,displacement):
        self.cppparticledata.setDisplacement(displacement)
    
    def get_displacement(self):
        return self.cppparticledata.getDisplacement()

    def set_diameters(self,diameter):
        self.cppparticledata.setDiameter(diameter.astype('float64'))
    
    def get_diameters(self):
        return self.cppparticledata.getDiameter()
    
    def set_mass(self,mass):
        self.cppparticledata.setMass(mass.astype('float64'))
    
    def get_mass(self):
        return self.cppparticledata.getMass()
    
    def set_position(self,position):
        self.cppparticledata.setAtomPosition(position.astype('float64'))
    
    def get_position(self):
        return self.cppparticledata.getAtomPosition()
    
    def set_velocity(self,velocity):
        self.cppparticledata.setAtomVelocity(velocity.astype('float64'))
    
    def get_velocity(self):
        return self.cppparticledata.getAtomVelocity()
    
    def get_neighbors(self, point, radius):
        return self.cppparticledata.getNeighbors(point,radius)
    def get_neighborsdistance(self, tag, radius):
        return self.cppparticledata.getNeighborsDistance(tag,radius)
