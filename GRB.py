import math

import numpy as np
import scipy as sp
from scipy.integrate import quad
from astropy.cosmology import Planck15


class GRB():
    "This is a class for a GRB. It accept redshift, spectal index and other information to describe a GRB"
    def __init__(self, z = 0, alpha = False, beta = False, Ep = 200, Fen = False, Fpho = False, SpecType = False, rangeup = False, rangedown = False):
        self.z = z
        self.alpha = self.alpha
        self.beta = self.beta
        self.Ep = self.Ep
        self.Fen = Fen
        self.Fpho = Fpho
        self.SpecType = SpecType
        self.rangeup = rangeup
        self.rangedown = rangedown


    @property
    def z(self):
        return self.__z
    @z.setter
    def z(self, z):
        assert (z >= 0), "The redshift must larger or equal to zero"
        self.__z = z

    @property
    def alpha(self):
        return self.__alpha
    @alpha.setter
    def alpha(self, alpha):
        assert alpha, "Please give alpha or alpha can not be zero" 
        self.__alpha = alpha

    @property
    def Ep(self):
        return self.__Ep
    @Ep.setter
    def Ep(self, Ep):

    @property
    def SpecType(self):
        return self.__SpecType
    @SpecType.setter
    def SpecType(self, spectype):
        if spectype:
            assert (spectype in ['Band', 'ExpCutoff']), "This SpecType are not known."
            self.__SpecType = spectype
        else:
            if self.beta :
                self.__SpecType = 'Band'
            else:
                self.__SpecType = "ExpCutoff"

    @property
    def rangedwon(self):
        return self.__rangedown
    @rangedown.setter
    def rangedown(self, rangedown):
        assert (self.rangeup > rangedown), "The up range must be larger than down limit."
        self.__rangedown = rangedown
        



    def __Band(self, E, A):
        Eb = (self.alpha - self.beta) * self.Ep / (2 + self.alpha)
        if E < Eb:
            return A * (E / 100) ** self.alpha * math.exp(- (2 + self.alpha) / Ep)
        else:
            return A * (E / 100) ** self.beta * math.exp((Eb / 100) ** (self.alpha - self.beta)) ** (self.beta - self.alpha) 

    def __ExpCutoff(self, E, A):
        return A * (E / 100) ** self.alpha ** self.exp(- (2 + self.alpha) * E / self.Ep)

    def Kcorrection(self, Eup, Edown):
        if self.SpecType == "Band":
            functemp = lambda E: E * self.__Band(self, E, 100)
        else self.SpecType == "ExpCutoff":
            functemp = lambda E: E * self.__SpecType(self, E, 100)
        downtemp = quad(functemp, self.rangedown, self.rangeup)
        uptemp = quad(functemp, Edown / (1 + self.z), Eup / (1 + self.z))
        self.Kcorr = uptemp / downtemp
        return self.Kcorr

    def L(self, cosmo=Planck15):
        if self.Fen:
            

