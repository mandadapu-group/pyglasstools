#implements the Simulation Box Class we've defined on the C++ side
from pyglasstools import _pyglasstools
import pyglasstools

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
