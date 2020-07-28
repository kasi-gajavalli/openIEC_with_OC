#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 15:06:08 2020

@author: kg245220
"""
import os
import sys
sys.path.append('/home/kg245220/code/pyOC')
sys.path.append('/home/kg245220/code/openiec_with_OC/openiec/property')
sys.path.append('/home/kg245220/code/openiec_with_OC/openiec/model')
sys.path.append('/home/kg245220/code/openiec_with_OC/openiec/calculate')
import pyOC
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
from pyOC import GridMinimizerStatus as gmStat
from molarVolume import MolarVolume_OC
from sigmacoint import SigmaCoherentInterface_OC
from minimize import SearchEquilibrium, ComputeEquilibrium
import numpy as np
from xarray import Dataset

pathname='/home/kg245220/code/OpenCalphad/cea_api_tests/Java/data/'
# Given temperature.
T= 3000
#Given pressure
P=1E5
# Render thermodynamic database.
db = 'feouzr.tdb'
 # Define components in the interface.
comps = ('U', 'O', 'ZR')
# Given initial alloy composition. x0 corresponds to the mole fractions of Al and Cr.
x0 = {
	'U' : 0.343298,
	'O' : 0.414924,
	'ZR': 0.241778
} 


limit=[0, 1.0]
dx=0.01
# Two phases separated by the interface : format phasename:(0,1).
phasenames = ('LIQUID','C1_FCC')  

# mass density laws (from Barrachin2004)
coriumMassDensityLaws = {
	'U1'   : lambda T: 17270.0-1.358*(T-1408),
	'ZR1'  : lambda T: 6844.51-0.609898*T+2.05008E-4*T**2-4.47829E-8*T**3+3.26469E-12*T**4,
	'O2U1' : lambda T: 8860.0-9.285E-1*(T-3120),
	'O2ZR1': lambda T: 5150-0.445*(T-2983),
	'O1'   : lambda T: 1.141 # set to meaningless value but ok as, no 'free' oxygen in the considered mixtures
}

##set verbosity
setverb=False

## tdb filepath
tdbFile=os.environ.get('OCDATA')+'/feouzr.tdb'
## reading tdb
elems=('O', 'U', 'ZR')
oc.readtdb(tdbFile,elems)

# set pressure
oc.setPressure(1E5)

# set temperature
oc.setTemperature(2800)
## suspend all phases except the liquid one
oc.setPhasesStatus(('* ',),phStat.Suspended)
oc.setPhasesStatus(('LIQUID',),phStat.Entered)

oc.setElementMolarAmounts(x0)

# calculate equilibrium
oc.changeEquilibriumRecord('eq2')
oc.calculateEquilibrium(gmStat.Off)

# retrieving mu data
phasesAtEquilibrium=oc.getPhasesAtEquilibrium()
phaseConstituentComposition = phasesAtEquilibrium.getPhaseConstituentComposition() 
mueq_inf=oc.getChemicalPotentials()


# calculate equilibrium with the grid-minimizer (equilibrium record = default equilibrium)
oc.changeEquilibriumRecord()
oc.calculateEquilibrium()

# retrieving mu data
phasesAtEquilibrium=oc.getPhasesAtEquilibrium()
phaseConstituentComposition = phasesAtEquilibrium.getPhaseConstituentComposition() 
mueq=oc.getChemicalPotentials()

# Use of molarVolume class
epsilon=0.1
molar_vol = MolarVolume_OC()
partialMolarVolumes,exactVolume,approxVolume=molar_vol.calculatePartialMolarVolume(T,x0,coriumMassDensityLaws,epsilon)

print(partialMolarVolumes,'\n')


sigma_model = SigmaCoherentInterface_OC(mueq_inf, mueq, partialMolarVolumes)

components = [each for each in elems]
cum = int(len(components))
print(
    "\n******************************************************************************\nOpenIEC is looking for interfacial equilibirium coposition.\nFor more information visit http://....../openiec."
)
x_s = SearchEquilibrium(sigma_model.objective, [limit] * cum, [dx] * cum)
x_c = ComputeEquilibrium(sigma_model.objective, x_s["x"])
print(
    "******************************************************************************\n\n"
)

print(x_c)

#sigma = sigma_model.infenergy(x_c)
#
#xx0 = [1.0 - sum(list(x0))] + list(x0)
#xx_c = [1.0 - sum(list(x_c))] + list(x_c)
#sigmapartial = list(np.array(sigma).flatten())
#sigmaavg = np.average([each for each in sigma])
#
#
#res = Dataset(
#        {
#            "Components": elems,
#            "Temperature": T,
#            "Initial_Alloy_Composition": ("Components", xx0),
#            "Interfacial_Composition": ("Components", xx_c),
#            "Partial_Interfacial_Energy": ("Components", sigmapartial),
#            "Interfacial_Energy": sigmaavg,
#        }
#    )
#
#print(res,'\n')






















