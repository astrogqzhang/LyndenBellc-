import math

import numpy as np
import scipy as sp
import scipy.interpolate as itp
from astropy.cosmology import Planck15 
 
class star():
    "This class is the base class for Limit and LyndenBell. It give the redshift and luminosity for an object"
    def __init__(self, z, L):
        self.z = z
        self.L = L
    
    @property
    def z(self):
        return self.__z
    @z.setter
    def z(self, z):
        assert (z > 0).all(), "redshift must be nonzero and non-negative"
        self.__z = z

    @property
    def L(self):
        return self.__L
    @L.setter
    def L(self, L):
        assert (L > 0).all(), "Luminosity must be nonzero and non-negative"


class Limit(star):
    "This class is prepared for LyndenBell. It return the flux limit of L and the limit of z."
    def __init__(self, z, L, xlim = 0, ylim = 0, Flim=0, cosmo=Planck15):
        super().__init__(z, L)
        self.xlim = xlim
        self.ylim = ylim
        self.Flim = Flim 
        self.cosmo = Planck15

    @property
    def xlim(self):
        return self.__xlim

    @xlim.setter
    def xlim(self):


class LyndenBell():
    """A class with redshift and luminosity. It want to compute the luminosity function and redshift for some objects with Lynden-Bell c- method.
    It also consider the luminosity evolution (1+z)^k."""
    def __init__(self, z, L, k = 0):
        super().__init__(z, L)
        self.k = k
        self.L0 = self.L / (1 + self.z) ** self.k

    @property
    def k(self):
        return self.__k

    @k.setter
    def k(self, k):
        self.__k = k



    
    
