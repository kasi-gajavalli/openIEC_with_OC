from coherentenergy import CoherentGibbsEnergy
from molarvolume import MolarVolume, InterficialMolarVolume
from openiec.calculate.calcsigma import SigmaCoherent
from openiec.calculate.minimize import SearchEquilibrium, ComputeEquilibrium
from pycalphad import Database
from openiec.model.sigmacoint import SigmaCoherentInterface
'''
INPUT
'''
# Given temperature.
T = 800.00
# Given initial alloy composition. x0 is the mole fraction of Al.
x0 = [0.2]  
# Render thermodynamic database.
db = Database("NiAlHuang1999.tdb")
# Define components in the interface.
comps = ["NI", "AL", "VA"]
# Two phases separated by the interface.
phasenames = ["FCC_A1", "GAMMA_PRIME"]
# Molar volumes of pure components to construct corresponding molar volume database.
# Molar volume of Ni.
vni = "6.718*10.0**(-6.0) + (2.936*10.0**(-5)*10.0**(-6.0))*T**1.355" 
# Molar volume of Al.
val = "10.269*10.0**(-6.0) + (3.860*10.0**(-5)*10.0**(-6.0))*T**1.491"
purevms = [[vni, val], ]*2

# A composition range for searching initial interfacial equilirium composition.
limit = [0.0001, 0.3]
# The composition step for searching initial interfacial equilirium composition.
dx = 0.1

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



print('**********************************************************************************', '\n')
'''------------------------------------------------------------
MOLAR VOLUME CALCUATIONS
------------------------------------------------------------'''
alphaFunc = MolarVolume(db, phasenames[0], comps, purevms[0])

betaFunc = MolarVolume(db, phasenames[1], comps, purevms[1])

vm = InterficialMolarVolume(alphaFunc, betaFunc)

print('\n')
print('The molar volume of bulk phase (alpha): ', vm[0]([0.9], 100.0), '\n')
print('\n')
print('The molar volume of bulk phase (beta): ', vm[1]([0.9], 100.0), '\n')

'''------------------------------------------------------------
SIGMA COHERENT CALCUATIONS
------------------------------------------------------------'''

sigma =  SigmaCoherent(T, x0, db, comps, phasenames, purevms, intervms=[], limit=[0, 1.0], dx=0.01)

print(sigma)




