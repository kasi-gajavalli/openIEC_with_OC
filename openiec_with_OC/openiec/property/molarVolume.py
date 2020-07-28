#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 10:49:07 2020

@author: kg245220
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import pyOC
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
from pyOC import GridMinimizerStatus as gmStat


import sympy as sy
import pycalphad.variables as V
from pycalphad import Database, Model
from sympy import lambdify, symbols, sympify, diff
import itertools
from functools import reduce


class MolarVolume_OC(object):
        
    def convertConstituentMolarToMassFractions(self,constituentMolarFractions,constituentsDescription):
    	constituentMassFractions=constituentMolarFractions.copy()
    	tot=0.0
    	for constituent in constituentMassFractions:
    		constituentMassFractions[constituent] *= constituentsDescription[constituent]['mass']
    		tot += constituentMassFractions[constituent]
    	fac = 1.0/tot
    	for constituent in constituentMassFractions:
    		constituentMassFractions[constituent] *= fac
    	return constituentMassFractions
    	
    
    ## function to calculate molar volume from constituent 
    def calculateVolume(self,temperature,phaseConstituentComposition,constituentsDescription,constituentMassDensityLaws):
    		print(phaseConstituentComposition.keys())
    		if (len(phaseConstituentComposition) != 1):
    			print('error: not a single phase (%s) at equilibrium for molar volume calculation' % list(phaseConstituentComposition.keys()))
#    			sys.exit()
    		constituentMolarFractions=phaseConstituentComposition[list(phaseConstituentComposition.keys())[0]]
    		# mass fractions from molar fractions
    		constituentMassFractions=self.convertConstituentMolarToMassFractions(constituentMolarFractions,constituentsDescription)
    		# ideal mixing law to evaluate mixture density
    		density=0.0
    		for constituent, massFraction in constituentMassFractions.items():
    			density += massFraction/constituentMassDensityLaws[constituent](temperature)
    		density=1.0/density
    		# total mass (mass is 'B' in OC)
    		mass=oc.getScalarResult('B')
    		return mass/density
    		
    		
    ## evaluate partial molar volumes by an approximation of the first order volume derivative by a second-order finite difference formula
    def calculatePartialMolarVolume(self,temperature,elementMolarAmounts,constituentMassDensityLaws,epsilon):
    	# set temperature
    	oc.setTemperature(temperature)
    	# evaluate (elementwise) partial molar volume (approximation of first order volume derivative by a second-order finite difference formula)
    	partialMolarVolumes={}
    	for element in elementMolarAmounts:
    		# evaluate volume for n[element]+epsilone
    		modifiedElementMolarAmounts=elementMolarAmounts.copy()
    		modifiedElementMolarAmounts[element] += epsilon
    		oc.setElementMolarAmounts(modifiedElementMolarAmounts)
    		oc.calculateEquilibrium(gmStat.Off)
    		volumePlus = self.calculateVolume(temperature,oc.getPhasesAtEquilibrium().getPhaseConstituentComposition(),oc.getConstituentsDescription(),constituentMassDensityLaws)
    		# evaluate volume for n[element]-epsilone
    		modifiedElementMolarAmounts[element] -= 2.0*epsilon
    		oc.setElementMolarAmounts(modifiedElementMolarAmounts)
    		oc.calculateEquilibrium(gmStat.Off)
    		volumeMinus = self.calculateVolume(temperature,oc.getPhasesAtEquilibrium().getPhaseConstituentComposition(),oc.getConstituentsDescription(),constituentMassDensityLaws)
    		partialMolarVolumes[element]=(volumePlus-volumeMinus)/(2.0*epsilon)
    	# verify that V = sum_i n_i V_i
    	oc.setElementMolarAmounts(elementMolarAmounts)
    	oc.calculateEquilibrium(gmStat.Off)
    	exactVolume = self.calculateVolume(temperature,oc.getPhasesAtEquilibrium().getPhaseConstituentComposition(),oc.getConstituentsDescription(),constituentMassDensityLaws)
    	approxVolume = 0.0
    	for element, molarAmount in oc.getPhasesAtEquilibrium().getPhaseElementComposition()[list(oc.getPhasesAtEquilibrium().getPhaseElementComposition())[0]].items():
    		approxVolume+=molarAmount*partialMolarVolumes[element]
    	return partialMolarVolumes,exactVolume,approxVolume
    
class MolarVolume(object):
    """
    Construct molar volumes of bulk phases.

    Parameters
    -----------
    db : Database
        Database containing the relevant parameters.
    comps : list
        Names of components to consider in the calculation.
    phasename : str
        Names of the phase to consider in the calculation.    
    purevm: list 
        The molar volume of the components.
    vm: a sympy expression
        The molar volume of the bulk phase.
    """

    def __init__(self, db, phasename, comps, purevm, intervm=[]):
        self.xs = [V.X(each) for each in comps if each != "VA"]
        self.vars_xs = [
            (self.xs[0], 1.0 - sum([self.xs[i] for i in range(1, len(self.xs))]))
        ]
        self.xxs = [self.xs[i] for i in range(1, len(self.xs))]
        self.vm = reduce(
            lambda x, y: x + y, [x * sympify(v) for x, v in zip(self.xs, purevm)]
        ).subs({"T": V.T})


def InterficialMolarVolume(alphavm, betavm):
    """
    Construct the partial molar volume of the interface.

    Parameters
    -----------
    alphavm: a sympy expression
        The molar volume of a bulk phase.
    betavm: a sympy expression
        The molar volume of another bulk phase.
    vm: a sympy expression
        The molar volume of the interfacial layer.
    vmis: list
        The partial molar volumes of components in the interface.
    """
    xs = alphavm.xs
    vm = 0.5 * (alphavm.vm + betavm.vm)
    dvmdxs = [diff(vm, x) for x in xs]

    sumvmi = reduce(lambda x, y: x + y, [x * dvmdx for x, dvmdx in zip(xs, dvmdxs)])

    vmis = [vm + dvmdx - sumvmi for dvmdx in dvmdxs]

    return [lambdify((alphavm.xxs, V.T), vmi, "numpy", dummify=True) for vmi in vmis]


if __name__ == "__main__":
    db = Database("NiAl.tdb")
    alphavm = MolarVolume(db, "FCC_A1", ["Ni", "Al"], ["1.0*T", "2.0*T"])
    betavm = MolarVolume(db, "LIQUID", ["Ni", "Al"], ["1.0*T", "2.0*T"])
    vm = InterficialMolarVolume(alphavm, betavm)

    print(vm[0]([0.9], 100.0), vm[1]([0.9], 100.0))

