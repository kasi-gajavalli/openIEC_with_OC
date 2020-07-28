#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 11:01:22 2020

@author: kg245220
"""

import os
import sys
sys.path.append('/home/kg245220/code/pyOC')
sys.path.append('/home/kg245220/code/openiec_with_OC/openiec/property')
import numpy as np
import matplotlib.pyplot as plt
import pyOC
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
#from pyOC import GridMinimizerStatus as gmStat
from molarVolume import MolarVolume_OC

oc.setVerbosity(False)
## tdb filepath
tdbFile=os.environ.get('OCDATA')+'/feouzr.tdb'
## reading tdb
elems=('O', 'U', 'ZR')
oc.readtdb(tdbFile,elems)
## suspend all phases except the liquid one
oc.setPhasesStatus(('* ',),phStat.Suspended)
oc.setPhasesStatus(('LIQUID',),phStat.Entered)
## set pressure
oc.setPressure(1E5)

# mass density laws (from Barrachin2004)
coriumMassDensityLaws = {
	'U1'   : lambda T: 17270.0-1.358*(T-1408),
	'ZR1'  : lambda T: 6844.51-0.609898*T+2.05008E-4*T**2-4.47829E-8*T**3+3.26469E-12*T**4,
	'O2U1' : lambda T: 8860.0-9.285E-1*(T-3120),
	'O2ZR1': lambda T: 5150-0.445*(T-2983),
	'O1'   : lambda T: 1.141 # set to meaningless value but ok as, no 'free' oxygen in the considered mixtures
}

# temperature and composition for which partial molar volumes are to be evaluated
temperature=3000
elementMolarAmounts = {
	'U' : 0.343298,
	'O' : 0.414924,
	'ZR': 0.241778
}

# testing for different values of the epsilon parameter that is used in the finite difference formula)
epsilons=np.array([1E-1,1E-2,1E-3,1E-4,1E-5,1E-6])
volumeErrors=np.empty(epsilons.size)
partialMolarVolumesU=np.empty(epsilons.size)
partialMolarVolumesZR=np.empty(epsilons.size)
partialMolarVolumesO=np.empty(epsilons.size)
for i in range(epsilons.size):
    epsilon=epsilons[i]
    molar_vol = MolarVolume_OC()
    partialMolarVolumes,exactVolume,approxVolume=molar_vol.calculatePartialMolarVolume(temperature,elementMolarAmounts,coriumMassDensityLaws,epsilon)
    volumeErrors[i]=(approxVolume/exactVolume-1.0)*100.0
    partialMolarVolumesU[i]=partialMolarVolumes['U']
    partialMolarVolumesZR[i]=partialMolarVolumes['ZR']
    partialMolarVolumesO[i]=partialMolarVolumes['O']
    print("partial molar volumes: ",partialMolarVolumes)
    print('volume error = {0:3.2e}%  '.format(volumeErrors[i]))
	
plt.rcParams["figure.figsize"] = (12,7)
fig,ax=plt.subplots(2,2,constrained_layout=True)
cax = ax[0,0]
csf = cax.loglog(epsilons,np.abs(volumeErrors), color='red',linestyle='dashed',marker='o')
cax.set_xlabel("$\\varepsilon$",fontsize=12)
cax.set_ylabel("|volume error| (%)",fontsize=12)
cax.grid(True)
cax = ax[0,1]
csf = cax.loglog(epsilons,partialMolarVolumesU, color='blue',linestyle='dashed',marker='x')
cax.set_xlabel("$\\varepsilon$",fontsize=12)
cax.set_ylabel("$V_U$",fontsize=12)
cax.grid(True)
cax = ax[1,0]
csf = cax.loglog(epsilons,partialMolarVolumesZR, color='blue',linestyle='dashed',marker='x')
cax.set_xlabel("$\\varepsilon$",fontsize=12)
cax.set_ylabel("$V_{Zr}$",fontsize=12)
cax.grid(True)
cax = ax[1,1]
csf = cax.loglog(epsilons,partialMolarVolumesO, color='blue',linestyle='dashed',marker='x')
cax.set_xlabel("$\\varepsilon$",fontsize=12)
cax.set_ylabel("$V_O$",fontsize=12)
cax.grid(True)
plt.show()