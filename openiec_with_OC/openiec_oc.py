"""According to given conditions and input parameters, calculate corresponding interfacial energies.
"""
import os
import pyOC
import numpy as np
import matplotlib.pyplot as plt
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
from pyOC import GridMinimizerStatus as gmStat
import math

######
# An example of Interfacil energy calculation
######

# oc setup
## setting verbosity (True or False - default), if set to yes, in particular, when getters are called the returned values are displayed in a comprehensive way
oc.setVerbosity(False)
## tdb filepath
tdbFile=os.environ.get('OCDATA')+'/feouzr.tdb'
## reading tdb
elems=('O', 'U', 'ZR')
oc.readtdb(tdbFile,elems)
## suspend all phases except the liquid one
oc.setPhasesStatus(('*',),phStat.Suspended)
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

'''------------------------------------------------------------
COHERENT GIBBS ENERGY CALCUATIONS
------------------------------------------------------------'''

coherentEnergy_obj = CoherentGibbsEnergy(T, db, comps, phasenames)

calculated_phase_equilibrium = coherentEnergy_obj.eqfunc(x0)
name_of_phase_in_equilibrium = coherentEnergy_obj.phase(x0)
phasevertex = coherentEnergy_obj.phasevertex(x0)
gm_of_bulk_phases = coherentEnergy_obj.Gibbsenergy(x0)
mu_of_components_in_bulk_phases = coherentEnergy_obj.chemicalpotential(x0)
phase_fractions = coherentEnergy_obj.phasefraction(x0)
mole_fractions = coherentEnergy_obj.molefraction(x0)
componentindex = coherentEnergy_obj.componentindex(x0)
site_fractions = coherentEnergy_obj.sitefraction(x0)
#

print('The calculated equilibrium is: ', calculated_phase_equilibrium, '\n')
print('The string name of the phase in equilibrium with the given conditions is: ', name_of_phase_in_equilibrium, '\n')
print('Phasevertex is the index of the phase in equilibrium: ', phasevertex, '\n')
print('Molar Gibbs energy of bulk phases: ', gm_of_bulk_phases, '\n')
print('Chemical potentials of components in bulk phases: ', mu_of_components_in_bulk_phases, '\n')
print('Phase fractions of phases in equilibrium: ', phase_fractions, '\n')
print('Mole fractions of components for phases in equilibrium: ', mole_fractions, '\n')
print('The componentindex is the index of the component in equilibrium: ', componentindex, '\n')
print('The Site fractions of components for phases in equilibrium: ', site_fractions, '\n')

"""
##Calculate the phase equilibrium, Gibbs energy and Chemical potentials of the associated components in the bulk phases
class CoherentGibbsEnergy(object):
   
    ##Calculate equilibrium of the bulk phases

    def Gm_bulkphase(temperature,elementMolarAmounts):
	# set temperature
    	oc.setTemperature(temperature)
        setElementMolarAmounts(elementMolarAmounts)
        bulkphase={}	
        bulkphase = oc.calculateEquilibrium(gmStat.On)
        c.getScalarResult('G')
        se_description = oc.getPhasesAtEquilibrium().getPhaseConstituentComposition()
        micalpotential=oc.getComponentAssociatedResult('MU')
"""
