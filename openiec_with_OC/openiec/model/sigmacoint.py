"""
Contruct the partial interfacial energies for different components on the solid/liquid interface and the object function to resolving the interfacial equilibrium condition.
"""

from openiec.utils.makemultigrid import makemultigrid
from scipy.optimize import minimize
from math import fabs
import numpy as np
from functools import reduce

class SigmaCoherentInterface(object):
    """
    Contruct the interfacial energy of the coherent interface and their partial quantaties.

    Parameters
    -----------
    alphafuncs: list of functions
        The phase with a single sublattice.
    betafuncs: tuple of functions
        The phase with two sublattices.
    mueq: list
        The chemical potentials of the equilibrium state at the given compositions.
    vmis: list
        The partial molar volumes of components.
    """

    def __init__(self, alphafuncs, betafuncs, mueq, vmis, comps):
        self.alphafuncs = alphafuncs
        self.betafuncs = betafuncs
        self.mueq = mueq
        self.vmis = vmis
        self.Nav = 6.02 * 10.0 ** (23.0)
        self.comps = comps

    def infenergy(self, x):
        """
        Compute coherent partial interfacial energies of various components.
        Parameters
        ----------
        x: list
            Interfacial composition.
        """
        x1 = list(x)
        
        x1 = {self.comps[0]:x[0], self.comps[1]:x[1], self.comps[2]:1-x[0]-x[1]}   # ASK ????????????????????????????????????
        
        sigma = [
            2.48
            * (
                0.5 * (self.alphafuncs(x1)[i] + self.betafuncs(x1)[i])
                - self.mueq[i]
            )
            * ((self.vmis[i](x) ** (-2.0 / 3.0)) * (self.Nav ** (-1.0 / 3.0)))
            for i in range(len(x) + 1)
        ]
        return sigma

    def objective(self, x):
        """
        Compute the absolute value of the differences between the partial interfacial energy.
        
        Parameters
        ----------
        x: list
            Interfacial composition.
        """
        v = 0.0
        s = self.infenergy(x)

        for i in range(len(s) - 1):
            for j in range(i + 1, len(s)):
                v += 10 * fabs(s[i] - s[j])
        return v