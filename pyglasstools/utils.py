from pyglasstools import _pyglasstools
import numpy as np
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
            self.SimBox = _pyglasstools.SimBox(L,origin,ndim);
        else:
            boxsize = np.array([Lx,Ly,Lz])#.astype('float64')
            self.SimBox = _pyglasstools.SimBox(boxsize,origin,ndim);
    
    #Redefine attributes so that it directly access SimBox C++ class 
    #attributes
    def __getattr__(self,attr):
            orig_attr = self.SimBox.__getattribute__(attr)
            if callable(orig_attr):
                def hooked(*args, **kwargs):
                    self.pre()
                    result = orig_attr(*args, **kwargs)
                    # prevent SimBox from becoming unwrapped
                    if result == self.SimBox:
                        return self
                    self.post()
                    return result
                return hooked
            else:
                return orig_attr
    
    ## \internal
    # \brief Get a C++ boxdim
    def _getSimBox(self):
        return self.SimBox

def read_gsd(filename,frame_num = 0):
    traj_gsd = gsd.hoomd.open(name=filename,mode='rb')
    return data_gsd(traj_gsd,frame_num)

class data_gsd(object):
    R""" Object to store particle data based on GSD format

    """
    def __init__(self, traj_gsd,frame_num):
        #Store trajectory data
        self.traj = traj_gsd; 
        #Store simulation box
        self.simbox = simbox(   Lx=self.traj[frame_num].configuration.box[0],
                                Ly=self.traj[frame_num].configuration.box[1],
                                Lz=self.traj[frame_num].configuration.box[2],
                                ndim=2)
        # create the c++ mirror class
        self.particledata = _pyglasstools.ParticleSystem(   self.simbox._getSimBox(), 
                                                            len(self.traj[frame_num].particles.diameter),
                                                            self.traj[frame_num].particles.diameter,
                                                            self.traj[frame_num].particles.mass,
                                                            self.traj[frame_num].particles.position,
                                                            self.traj[frame_num].particles.velocity 
                                                            );
    def update(self,frame_num):
        #Store simulation box
        self.simbox = simbox(   Lx=self.traj[frame_num].configuration.box[0],
                                Ly=self.traj[frame_num].configuration.box[1],
                                Lz=self.traj[frame_num].configuration.box[2],
                                ndim=2)
        
        # create the c++ mirror class
        self.particledata = _pyglasstools.ParticleSystem(   self.simbox._getSimBox(), 
                                                            len(self.traj[frame_num].particles.diameter),
                                                            self.traj[frame_num].particles.diameter,
                                                            self.traj[frame_num].particles.mass,
                                                            self.traj[frame_num].particles.position,
                                                            self.traj[frame_num].particles.velocity 
                                                            );
    def _getParticleSystem(self):
        return self.particledata
    
    def set_diameters(self,diameter):
        self.particledata.setDiameter(diameter.astype('float64'))
    
    def get_diameters(self):
        return self.particledata.getDiameter()
    
    def set_mass(self,mass):
        self.particledata.setMass(mass.astype('float64'))
    
    def get_mass(self):
        return self.particledata.getMass()
    
    def set_position(self,position):
        self.particledata.setAtomPosition(position.astype('float64'))
    
    def get_position(self):
        return self.particledata.getAtomPosition()
    
    def set_velocity(self,velocity):
        self.particledata.setAtomVelocity(velocity.astype('float64'))
    
    def get_velocity(self):
        return self.particledata.getAtomVelocity()
    
    def get_neighbors(self, point, radius):
        return self.particledata.getNeighbors(point,radius)
