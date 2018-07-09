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
        self.__L = L


class Limit(star):
    "This class is prepared for LyndenBell. It return the flux limit of L and the limit of z."
    def __init__(self, z, L, zlim = 0, Llim = 0, Flim=0, k = 0, cosmo=Planck15):
        super().__init__(z, L)
        self.zlim = xlim
        self.Llim = ylim
        self.Flim = Flim 
        self.k = k
        self.cosmo = cosmos 
        self.L0 = self.L / (1 + self.z) ** self.k
        if self.Flim != 0:
            __x = x.linspace(0, 10, 100)
            __dl = self.cosmos.luminosity_distance(__x).cgs.value
            __y = 4 * math.pi * __dl ** 2 * self.Flim / (1 + self.z) ** self.k


    @property
    def zlim(self):
        return self.__zlim
    @zlim.setter
    def zlim(self, zlim):
        if zlim == 0:
            assert Flim ！= 0， "zlim and Flim cannot both be zero"
            f = itp.interp1d(__y, __x)
            self.__zlim = f(self.L0)
        else: 
            f = itp.interp1d(self.Llim / (1 + self.z) ** self.k, self.z)
            self.__zlim = f(self.L0) 
    
    @property
    def Llim(self):
        return self.__Llim
    @Llim.setter
    def Llim(self, Llim):
        if Llim == 0:
            assert Flim != 0, "Llim and Flim cannot both be zero"
            self.__Llim = 4 * math.pi * self.cosmo.luminosity_distance(self.z).cgs.value ** 2 * self.Flim / (1 + self.z) ** self.k
        else
            self.__Llim = Llim / (1 + self.z) ** k



class LyndenBell():
    """A class with redshift and luminosity. It want to compute the luminosity function and redshift for some objects with Lynden-Bell c- method.
    It also consider the luminosity evolution (1+z)^k."""
    def __init__(self, z, L, lim, k = 0, cosmo = Planck15):
        super().__init__(z, L)
        self.k = k
        self.lim = lim
        self.cosmo = cosmo
        self.L0 = self.L / (1 + self.z) ** self.k

    def lyndenbellc(self):
        self.__Marr = []
        self.__Narr = []
        for i in range(len(self.x)):
            zijudge = self.z <= self.lim.zlim[i]
            Lijudge = self.L0 >= self.L0[i]
            judgearray = np.logical_and(zijudge, Lijudge)
            self.__Narr.append(judgearray.sum())
        for j in range(len(self.y)):
            zjjudge = self.z <= self.z[j]
            Ljjudge = self.L0 >= self.lim.Llim[j]
            judgearray = np.logical_and(zjjudge, Ljjudge)
            self.__Marr.append(judgearray.sum())
        zdis = []
        Ldis = []
        for i in range(len(x)):
            ztempfunc = 1
            Ltempfunc = 1
            ziarray = self.z < self.z[i]
            Liarray = self.L0 > self.L0[i]
            for j in range(len(self.z)):
                if ziarray[j] and self.__Marr[j]: ztempfunc *= (1 + 1 / self.__Marr[j])
            for j in range(len(self.L0)):
                if Liarray[j] and self.__Narr[j]: Ltempfunc *= (1 + 1 / self.__Narr[j])
            zdis.append(ztempfunc)
            Ldis.append(Ltempfunc)
        zdis = np.array(zdis)
        Ldis = np.array(Ldis)
        return zdis, Ldis

    def testindependence(self):
        __karray = []
        __tau = []
        for k in np.linspace(-5, 5, 2000):
            __L0 = self.L / (1 + self.z) ** k
            __Lim = Limit(self.z, self.L, self.lim.zlim, self.lim.Llim, self.lim.Flim, k, cosmo=self.cosmo)
            uptemp = 0
            downtemp = 0
            for i in range(len(self.x)):
                zijudge = self.z <= __lim.zlim[i]
                Lijudge = __L0 >= __L0[i]
                Ni = np.logical_and(zijudge, Lijudge)
                Rijudge = self.z <= self.z[i]
                Ri = np.logical_and(Rijudge, Lijudge)
                Ei = (1 + Ni) / 2
                Vi = (Ni ** 2 - 1) / 12
                uptemp += (Ri - Ei)
                downtemp += Vi 
            __tau.append(uptemp / math.sqrt(downtemp))
            __karray.append(k)
        __tau = np.array(__tau)
        __karray = np.array(__karray)
        print("The best fit of k is {}".format(__karray[np.abs(__tau.argmin())])) 
        print("The 1 sigma error is k is +{} -{}".format(__karray[np.abs(__tau -1).argmin()], __karray[np.abs(__tau + 1).argmin()]))
        self.k = __karray[np.abs(__tau.argmin())]



            

 
            
    


    


    
    
