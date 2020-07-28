import numpy as np
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
import matplotlib.pyplot as plt
from ver2_coherentenergy import CoherentGibbsEnergy_OC
from ver2_molarvolume import MolarVolume


def test():
    #print(os.getpwd()) 
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
   
    # Two phases separated by the interface : format phasename:(0,1).
    phasenames = ('LIQUID','C1_FCC')
    
    ##My addition
    
#    #Put 'True' to suspend all phases else 'False'
#    SuspendAllPhases=True
#    
#    #Phases to be entered
#    phase_entered = ('LIQUID')
#    ## suspend all phases except the liquid one
#    oc.setPhasesStatus(('* ',),phStat.Suspended)
#    oc.setPhasesStatus(('LIQUID',),phStat.Entered)
#    
    ## My addition completed
    
    ##set verbosity
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
    epsilons =np.array([1E-1,1E-2,1E-3,1E-4,1E-5,1E-6])
    
    partialMolarVolume_U=np.zeros(epsilons.size)
    partialMolarVolume_O=np.zeros(epsilons.size)
    partialMolarVolume_Z=np.zeros(epsilons.size)
    volumeErrors = np.zeros(epsilons.size)
    partialMolarVolumes={}
    
    coherentGibbsEnergy = CoherentGibbsEnergy_OC()
    coherentGibbsEnergy.readDatabase(pathname, db, comps)
    coherentGibbsEnergy.suspendPhases('* ', phasenames)
    
    molar_vol = MolarVolume()
    
    for i in range(epsilons.size):
     
        for elems in x0.keys():
                
            
            x0ep = x0.copy()
            
            x0ep[elems] += epsilons[i]
            coherentGibbsEnergy.calculateEquilibria('Off',setverb,T,P,x0ep)
            cd=coherentGibbsEnergy.getConstituentsDescription()
            pcc=coherentGibbsEnergy.getPhaseConstituentComposition()
            mass=coherentGibbsEnergy.getMass()
            
            x0ep[elems] -= 2.*epsilons[i]
            coherentGibbsEnergy.calculateEquilibria('Off',setverb,T,P,x0ep)
            cd_m=coherentGibbsEnergy.getConstituentsDescription()
            pcc_m=coherentGibbsEnergy.getPhaseConstituentComposition()
            mass_m=coherentGibbsEnergy.getMass()
            
            partial_molar_vol = molar_vol.calculatePartialMolarVolume(T,pcc,cd,pcc_m,cd_m,coriumMassDensityLaws, epsilons[i],mass,mass_m, phasenames[0])
            partialMolarVolumes[elems] = partial_molar_vol
            
        ###My addition##
        coherentGibbsEnergy.calculateEquilibria('Off',setverb,T,P,x0)
        
        cd=coherentGibbsEnergy.getConstituentsDescription()
        pcc=coherentGibbsEnergy.getPhaseConstituentComposition()
        pec= coherentGibbsEnergy.getPhaseElementComposition()
        mass=coherentGibbsEnergy.getMass()
        
        exactVolume, approxVolume = molar_vol.exact_and_approx_volume(T,pcc,cd,pec,coriumMassDensityLaws,mass,phasenames[0],partialMolarVolumes)
        
        volumeErrors[i]=(approxVolume/exactVolume-1.0)*100.0
        
        partialMolarVolume_U[i] = partialMolarVolumes['U']
        partialMolarVolume_O[i] = partialMolarVolumes['O']
        partialMolarVolume_Z[i] = partialMolarVolumes['ZR']
        
    return partialMolarVolume_U, partialMolarVolume_O, partialMolarVolume_Z, volumeErrors, epsilons


partialMolarVolumesU, partialMolarVolumesO, partialMolarVolumesZR, volumeErrors, epsilons =  test()

print(partialMolarVolumesU,partialMolarVolumesO,partialMolarVolumesZR)

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




