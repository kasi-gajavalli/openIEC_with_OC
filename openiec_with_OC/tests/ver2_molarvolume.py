"""
Construct the molar volume CALPHAD model.
"""
import os
import pyOC
import itertools
import numpy as np
import sympy as sy
import matplotlib.pyplot as plt
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
from pyOC import GridMinimizerStatus as gmStat
from sympy import lambdify, symbols, sympify, diff
from functools import reduce


class MolarVolume(object):

######
# Functions for partial molar volumes evaluation
######
## convert constituent molar fractions to mass fractions
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
    def calculateVolume(self,temperature,phaseConstituentComposition,constituentsDescription,constituentMassDensityLaws,mass, phase_name):
        if (len(phaseConstituentComposition) != 1):
            print('Warning: not a single phase (%s) at equilibrium for molar volume calculation' % list(phaseConstituentComposition.keys()))
#        for phases in phaseConstituentComposition.keys():
        constituentMolarFractions=phaseConstituentComposition[phase_name]
#        print('##################',phaseConstituentComposition[0],'#########################')
		# mass fractions from molar fractions
        constituentMassFractions=self.convertConstituentMolarToMassFractions(constituentMolarFractions,constituentsDescription)
		# ideal mixing law to evaluate mixture density
#            print('*******************',constituentMassFractions,'********************')
        density=0.0
        for constituent, massFraction in constituentMassFractions.items():
            density += massFraction/constituentMassDensityLaws[constituent](temperature)
        density=1.0/density
    
#		# total mass (mass is 'B' in OC)
        volume=mass/density
        return volume
    
    ## evaluate partial molar volumes by an approximation of the first order volume derivative by a second-order finite difference formula
    def calculatePartialMolarVolume(self,temperature,pcc,cd,pcc_m,cd_m,constituentMassDensityLaws, epsilon,mass,mass_m, phase_name):
        # evaluate (elementwise) partial molar volume (approximation of first order volume derivative by a second-order finite difference formula)
        
        # evaluate volume for n[element]+epsilone
        volumePlus = self.calculateVolume(temperature,pcc,cd,constituentMassDensityLaws,mass, phase_name)
        
        # evaluate volume for n[element]-epsilone
        volumeMinus = self.calculateVolume(temperature,pcc_m,cd_m,constituentMassDensityLaws,mass_m, phase_name)
        
        partialMolarVolumes = (volumePlus - volumeMinus)/(2.0*epsilon)
        
        return partialMolarVolumes
    
    def exact_and_approx_volume(self,T,pcc,cd,pec,constituentMassDensityLaws,mass,phase_name,partialMolarVolumes):
        
        exactVolume = self.calculateVolume(T,pcc,cd,constituentMassDensityLaws,mass,phase_name)
        approxVolume = 0.0
#        for element, molarAmount in oc.getPhasesAtEquilibrium().getPhaseElementComposition()[list(oc.getPhasesAtEquilibrium().getPhaseElementComposition())[0]].items():
        for element, molarAmount in pec[phase_name].items():
            approxVolume+=molarAmount*partialMolarVolumes[element]
            
        return exactVolume,approxVolume

#    epsilons=np.array([1E-1,1E-2,1E-3,1E-4,1E-5,1E-6])
#    volumeErrors=np.empty(epsilons.size)
#    partialMolarVolumesU=np.empty(epsilons.size)
#    partialMolarVolumesZR=np.empty(epsilons.size)
#    partialMolarVolumesO=np.empty(epsilons.size)
#    for i in range(epsilons.size):
#    	epsilon=epsilons[i]
#    	partialMolarVolumes,exactVolume,approxVolume=self.calculatePartialMolarVolume(temperature,elementMolarAmounts,constituentMassDensityLaws,epsilon)
#    	volumeErrors[i]=(approxVolume/exactVolume-1.0)*100.0
#    	partialMolarVolumesU[i]=partialMolarVolumes['U']
#    	partialMolarVolumesZR[i]=partialMolarVolumes['ZR']
#    	partialMolarVolumesO[i]=partialMolarVolumes['O']
#    	print("partial molar volumes: ",partialMolarVolumes)
#    	print('volume error = {0:3.2e}%  '.format(volumeErrors[i]))
