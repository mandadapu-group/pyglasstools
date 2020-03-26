#implements the Simulation Box Class we've defined on the C++ side
from pyglasstools import _pyglasstools
import pyglasstools
import numpy as np

class simbox(object):
    R""" Define box ndim.

    """
    def __init__(self, Lx=1.0, Ly=1.0, Lz=1.0, ndim=3, L=None, volume=None):
        if L is not None:
            Lx = L;
            Ly = L;
            Lz = L;

        if ndim == 2:
            Lz = 1.0;

        self.Lx = Lx;
        self.Ly = Ly;
        self.Lz = Lz;
        self.ndim = ndim;
        self.SimBox = _pyglasstools.SimBox(self.Lx, self.Ly, self.Lz, self.ndim);
        if volume is not None:
            self.set_volume(volume);

    ## \internal
    # \brief Get a C++ boxdim
    def _getSimBox(self):
        return self.SimBox

    def get_volume(self):
        R""" Get the box volume.

        Returns:
            The box volume (area in 2D).
        """
        return self.SimBox.getVolume();
    
    def get_dim(self):
        R""" Get the box volume.

        Returns:
            The box volume (area in 2D).
        """
        return self.SimBox.getDim();

    def apply_pbc(self,v):
        #Make sure that it's a float64 type
        return self.SimBox.applyPBC(v.astype('float64'))

    def __str__(self):
        return 'Simulation Box: Lx=' + str(self.Lx) + ' Ly=' + str(self.Ly) + ' Lz=' + str(self.Lz) + ' ndim=' + str(self.ndim);

class data(object):
    R""" Object to store particle data

    """
    def __init__(self, diameter, mass, position, velocity, simbox, mode=None):

        self.diameter = diameter
        self.mass = mass
        self.position = position
        self.velocity = velocity
        
        # create the c++ mirror class
        self.particledata = _pyglasstools.SystemData(len(diameter),self.diameter,self.mass,self.position,self.velocity, simbox._getSimBox());
    def _getSystemData(self):
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

class system(object):
    R""" Object to store system data, including particle and simulation box data

    """
    def __init__(self, simbox, sysdata):

        # create the c++ mirror class
        self.particlesystem = _pyglasstools.ParticleSystem(simbox._getSimBox(),sysdata._getSystemData());

    def get_neighborsid(self,position,radius):
        if (len(position) < 3):
            position = np.append(position,0)
        return self.particlesystem.getNeighborsID(position.astype('float64'),radius)
