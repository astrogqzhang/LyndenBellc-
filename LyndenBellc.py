import math
from itertools import repeat
from multiprocessing import Pool

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
    '''This class is prepared for LyndenBell. It return the flux limit of L and the limit of z. 
        Flim ----------- The flux limit. When zlim, Llim and Flim are both known, only conisder Flim.
        k -------------- The (1 + z) ^ k. This parameter is prepared for luminosity evolution.
        cosmo ---------- The cosmology model.
        Llim ----------- The luminosity limit. It is the minmum luminosity at z.
        zlim ----------- The redshift limit. It is the maximum redshift at which the GRB can be observed.'''
    def __init__(self, z, L, zlim = 0, Llim = 0, Flim=0, k = 0, cosmo=Planck15):
        super().__init__(z, L)
        self.Flim = Flim 
        self.k = k
        self.cosmo = cosmo 
        self.L0 = self.L / (1 + self.z) ** self.k
        self.__x = np.linspace(0, 10, 100)
        self.__dl = self.cosmo.luminosity_distance(self.__x).cgs.value
        self.__y = 4 * math.pi * self.__dl ** 2 * self.Flim / (1 + self.__x) ** self.k
        self.Llim = Llim
        self.zlim = zlim
        self.z_turnover = self.__x[self.__y.argmax()]
        if self.z_turnover == 0:
            self.z_turnover = self.z.max() + 1


    @property
    def zlim(self):
        return self.__zlim
    @zlim.setter
    def zlim(self, zlim):
        if self.Flim != 0:
            f = itp.interp1d(self.__y, self.__x, fill_value="extrapolate")
            self.__zlim = f(self.L0)
        else: 
            assert type(zlim) != type(0), "zlim and Flim can not both be zero"
            f = itp.interp1d(self.Llim / (1 + self.z) ** self.k, self.z, fill_value="extrapolate", kind='cubic')
            self.__zlim = f(self.L0)
            # If the luminosity larger than the maximum limit, then set this zlim to be maximum redshift plus 1.
            # This behavior is for avoid the error of interpolate.
        for i in range(len(self.L0)):
            if self.L0[i] > (self.Llim).max():
                self.__zlim[i] = self.z.max() + 1 
    
    @property
    def Llim(self):
        return self.__Llim
    @Llim.setter
    def Llim(self, Llim):
        if self.Flim != 0:
            self.__Llim = 4 * math.pi * self.cosmo.luminosity_distance(self.z).cgs.value ** 2 * self.Flim / (1 + self.z) ** self.k
        else:
            assert type(Llim) != type(0), "zlim and Flim can not both be zero"
            self.__Llim = Llim / (1 + self.z) ** self.k



class LyndenBell(star):
    """A class with redshift and luminosity. It want to compute the luminosity function and redshift for some objects with Lynden-Bell c- method.
    It also consider the luminosity evolution (1+z)^k."""
    def __init__(self, z, L, lim, k = 0, cosmo = Planck15):
        super().__init__(z, L)
        self.k = k
        self.lim = lim
        self.cosmo = cosmo
        self.L0 = self.L / (1 + self.z) ** self.k


    def lyndenbellc(self):
        self.lim = Limit(self.z, self.L, self.lim.zlim, self.lim.Llim, self.lim.Flim, self.k, cosmo=self.cosmo)
        idx = self.z < self.lim.z_turnover
        self.ztemp = self.z[idx]
        self.Ltemp = self.L[idx]
        zlim = self.lim.zlim[idx]
        Llim = self.lim.Llim[idx]
        self.L0 = self.Ltemp / (1 + self.ztemp) ** self.k
        self.__Marr = []
        self.__Narr = []
        for i in range(len(self.ztemp)):
            zijudge = self.ztemp <= zlim[i]
            Lijudge = self.L0 >= self.L0[i]
            judgearray = np.logical_and(zijudge, Lijudge)
            self.__Narr.append(judgearray.sum())
        for j in range(len(self.Ltemp)):
            zjjudge = self.ztemp <= self.ztemp[j]
            Ljjudge = self.L0 >= Llim[j]
            judgearray = np.logical_and(zjjudge, Ljjudge)
            self.__Marr.append(judgearray.sum())
        zdis = []
        Ldis = []
        for i in range(len(self.ztemp)):
            ztempfunc = 1
            Ltempfunc = 1
            ziarray = self.ztemp < self.ztemp[i]
            Liarray = self.L0 > self.L0[i]
            for j in range(len(self.ztemp)):
                if ziarray[j] and self.__Marr[j]: ztempfunc *= (1 + 1 / self.__Marr[j])
            for j in range(len(self.L0)):
                if Liarray[j] and self.__Narr[j]: Ltempfunc *= (1 + 1 / self.__Narr[j])
            zdis.append(ztempfunc)
            Ldis.append(Ltempfunc)
        zdis = np.array(zdis)
        Ldis = np.array(Ldis)
        return zdis, Ldis

    def testindependence(self, Weight='sqrtVar', processnumber = 4):
        """This function is written for test dependence between redshift and luminosity. It use tau static to 
        test dependence. It can receive a parameter as Width. This function only consider 'sqrtVar' and 'equal' as width."""
        Flag = Weight in ['sqrtVar', 'equal']
        assert Flag, "The width must be sqrtVar or equal"
        __tau = []
        # for k in np.linspace(-5, 5, 5000):
        #     tau = self.tau(k, Weight = Weight)
        #     __tau.append(tau)
        #     __karray.append(k)
        __karray = np.linspace(-5, 5, 5000)
        with Pool(processes = processnumber) as pool:
            __tau.append(pool.starmap(self.tau, zip(__karray, repeat(Weight))))
        __tau = np.array(__tau)
        print("The best fit of k is {}".format(__karray[np.abs(__tau).argmin()])) 
        print("The 1 sigma error is k is +1 sigma:{} -1 sigma:{}".format(__karray[np.abs(__tau -1).argmin()], __karray[np.abs(__tau + 1).argmin()]))
        self.k = __karray[np.abs(__tau).argmin()]
        return __karray, __tau
    
    def tau(self, k = 0, Weight = 'sqrtVar'):
        "This function consider (1+z)^k as the function to remove luminosity evolution."
        Lim = Limit(self.z, self.L, self.lim.zlim, self.lim.Llim, self.lim.Flim, k, cosmo = self.cosmo)
        idx = self.z <= Lim.z_turnover
        self.ztemp = self.z[idx]
        self.Ltemp = self.L[idx]
        zlim = Lim.zlim[idx]
        Llim = Lim.Llim[idx]
        L0 = self.Ltemp / (1 + self.ztemp) ** k
        Viarr = []
        Tiarr = []
        for i in range(len(self.ztemp)):
            zijudge = self.ztemp <= zlim[i]
            Lijudge = L0 >= L0[i]
            Ni = np.logical_and(zijudge, Lijudge).sum()
            if Ni == 1: continue # don't consider Ni == 1
            Rijudge = self.ztemp <= self.ztemp[i]
            Ri = np.logical_and(Rijudge, Lijudge).sum()
            Ei = (1 + Ni) / 2
            Vi = (Ni ** 2 - 1) / 12
            Ti = (Ri - Ei) / math.sqrt(Vi)
            Tiarr.append(Ti)
            Viarr.append(Vi)
        Viarr = np.array(Viarr)
        Tiarr = np.array(Tiarr)
        if Weight == 'sqrtVar':
            tautemp = (Tiarr * np.sqrt(Viarr)).sum() / math.sqrt(Viarr.sum()) 
        elif Weight == 'equal':
            tautemp = Tiarr.sum()
        return tautemp
