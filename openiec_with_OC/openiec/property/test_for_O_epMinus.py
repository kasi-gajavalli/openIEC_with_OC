#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 10:47:57 2020

@author: kg245220
"""

import sys
sys.path.append('/home/kg245220/code/pyOC')
import os
import pyOC
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
from pyOC import GridMinimizerStatus as gmStat
from coherentenergy import CoherentGibbsEnergy
from molarvolume import MolarVolume



pathname='/home/kg245220/code/OpenCalphad/cea_api_tests/Java/data/'
    # Given temperature.
T= 2773
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
   
# Two phases separated by the interface.
#phasenames = ('C1_FCC', 'LIQUID')
##set verbolsity
setverb=False



x0['O'] += 0.1
x0['O'] -= 2.0*0.1

# setting verbosity (True or False - default), if set to yes, in particular, when getters are called the returned values are displayed in a comprehensive way
oc.setVerbosity(True)

# tdb filepath
tdbFile=pathname+'feouzr.tdb'
print(tdbFile)
exit

# reading tdb
elems=('O', 'U', 'ZR')
oc.readtdb(tdbFile,elems)

# play with phase status
oc.setPhasesStatus(('C1_FCC',),phStat.Suspended)
phaseNames=('LIQUID',)
oc.setPhasesStatus(phaseNames,phStat.Entered)

# set pressure
oc.setPressure(1E5)

# set temperature
oc.setTemperature(2773)

'''------------------------------------------------------------
COHERENT GIBBS ENERGY CALCUATIONS
------------------------------------------------------------'''

oc.setElementMolarAmounts(x0)

# calculate equilibrium without the grid-minimizer (equilibrium record = eq2)
oc.calculateEquilibrium(gmStat.Off)

phasesAtEquilibrium=oc.getPhasesAtEquilibrium()

pcc = phaseConstituentComposition=phasesAtEquilibrium.getPhaseConstituentComposition()
mass=oc.getScalarResult('B')

#
#
#cge = CoherentGibbsEnergy(T, P, db, comps, x0, phasenames, setverb, pathname)
#cge.eqfunc('Off')

#cd=cge.getConstituentsDescription()
#pcc=cge.getPhaseConstituentComposition()
#mass=cge.getMass()

print(mass)