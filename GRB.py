import math

import numpy as np
import scipy as sp
from scipy.integrate import quad
from astropy.cosmology import Planck15, z_at_value
import astropy.units as u


class GRB():
    "This is a class for a GRB. It accept redshift, spectal index and other information to describe a GRB"
    def __init__(self, z = 0, L = 0, alpha = False, beta = False, Ep = 200, Fen = False, Fpho = False, SpecType = False, rangeup = False, rangedown = False):
        self.z = z
        self.L = L
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
        assert (Ep > 0), "The peak energy must larger than zero"
        self.__Ep = Ep

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
            return A * (E / 100) ** self.alpha * math.exp(- (2 + self.alpha) / self.Ep)
        else:
            return A * (E / 100) ** self.beta * math.exp((Eb / 100) ** (self.alpha - self.beta)) ** (self.beta - self.alpha) 

    def __ExpCutoff(self, E, A):
        return A * (E / 100) ** self.alpha ** math.exp(- (2 + self.alpha) * E / self.Ep)

    def Kcorrection(self, Eup, Edown, FLAG = False):
        if self.SpecType == "Band":
            functemp = lambda E: E * self.__Band(E, 100)
        elif self.SpecType == "ExpCutoff":
            functemp = lambda E: E * self.__SpecType(E, 100)
        if FLAG == 'photon':
            if self.SpecType == "Band":
                functemp2 = lambda E: self.__Band(E, 100)
            elif self.SpecType == "ExpCutoff":
                functemp2 = lambda E: self.__SpecType(E, 100)
            downtemp = quad(functemp2, self.rangedown, self.rangeup)
        else:
            downtemp = quad(functemp, self.rangedown, self.rangeup)
        uptemp = quad(functemp, Edown / (1 + self.z), Eup / (1 + self.z))
        self.Kcorr = uptemp / downtemp
        return self.Kcorr

    def obtainL(self, cosmo=Planck15, Fen = False, Fpho = False, ifKcorr = False):
        for i in [Fen, Fpho, self.Fen, self.Fpho]:
            if i:
                if i in [Fen, self.Fen]:
                    temp = i
                else:
                    temp = i * self.Kcorrection(self.rangeup, self.rangedown, FLAG='photon') 
        if not ifKcorr:
            try:
                temp = temp * self.Kcorr
            except:
                print("The K-correction is not adopted, Because K-correction isn't assigned.")
        return 4 * math.pi * cosmo.luminosity_distance(self.z).value.cgs ** 2 * temp
    
    def obtainz(self, cosmo=Planck15, Fen = False, Fpho = False, ifKcorr = False):
        for i in [Fen, Fpho, self.Fen, self.Fpho]:
            if i:
                if i in [Fen, self.Fen]:
                    temp = i
                else:
                    temp = i * self.Kcorrection(self.rangeup, self.rangedown, FLAG='photon') 
        if not ifKcorr:
            try:
                temp = temp * self.Kcorr
            except:
                print("The K-correction is not adopted, Because K-correction isn't assigned.")
        dl = math.sqrt(self.L / 4 / math.pi / temp)
        dl = dl * u.cm
        return z_at_value(cosmo.luminosity_distance, dl)
        
        
        
 


            

