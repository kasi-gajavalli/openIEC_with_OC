#from openiec.calculate.calcsigma import SigmaCoherent
import sys
sys.path.append('/home/kg245220/code/pyOC')
sys.path.append('/home/kg245220/code/openiec_with_OC/openiec/property')
import os
#import pyOC
#from pyOC import opencalphad as oc
#from pyOC import PhaseStatus as phStat
#from pyOC import GridMinimizerStatus as gmStat
from coherentenergy import CoherentGibbsEnergy
from molarvolume import MolarVolume



def test():
    #print(os.getpwd()) 
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
    phasenames = ('LIQUID','C1_FCC')
    ##set verbolsity
    setverb=False
         
#    cge = CoherentGibbsEnergy(T, P, db, comps, x0, phasenames,setverb,pathname)
#    cge.eqfunc('On')
#    cge.generate_phase_info()    
#    pec=cge.getPhaseElementFraction()
#    pcc=cge.getPhaseConstituentComposition()
#    cge.getPhaseSites()
#    cge.getChemicalPotential()
#    cd=cge.getConstituentsDescription()
#    ma=cge.getMolarAmounts()
    
#    print(ma)
    
 
    # mass density laws (from Barrachin2004)
    coriumMassDensityLaws = {
	'U1'   : lambda T: 17270.0-1.358*(T-1408),
	'ZR1'  : lambda T: 6844.51-0.609898*T+2.05008E-4*T**2-4.47829E-8*T**3+3.26469E-12*T**4,
	'O2U1' : lambda T: 8860.0-9.285E-1*(T-3120),
	'O2ZR1': lambda T: 5150-0.445*(T-2983),
	'O1'   : lambda T: 1.141 # set to meaningless value but ok as, no 'free' oxygen in the considered mixtures
}
    epsilon = 1E-1
    
    partialMolarVolumes={}
     
    for elems in x0.keys():
            
        
        x0ep = x0.copy()
        
        
        x0ep[elems] += epsilon
        cge = CoherentGibbsEnergy(T, P, db, comps, x0ep, phasenames, setverb, pathname)
        cge.eqfunc('Off',elems,False)
        cd=cge.getConstituentsDescription()
        pcc=cge.getPhaseConstituentComposition()
        mass=cge.getMass()
                
        
        x0ep[elems] -= 2.*epsilon
        cge_m = CoherentGibbsEnergy(T, P, db, comps, x0ep, phasenames,setverb,pathname)
        cge_m.eqfunc('Off',elems,False)
        cd_m=cge_m.getConstituentsDescription()
        pcc_m=cge_m.getPhaseConstituentComposition()
        mass_m=cge_m.getMass()
        
        molar_vol = MolarVolume()
        partial_molar_vol = molar_vol.calculatePartialMolarVolume(T,pcc,cd,pcc_m,cd_m,coriumMassDensityLaws, epsilon,mass,mass_m, phasenames[0])
        partialMolarVolumes[elems] = partial_molar_vol
    
    
    return partialMolarVolumes

if __name__ == "__main__":
    print(test())