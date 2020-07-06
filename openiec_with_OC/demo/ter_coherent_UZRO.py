#from openiec.calculate.calcsigma import SigmaCoherent
import sys
sys.path.append('/home/kg245220/code/pyOC')
import os
#import pyOC
#from pyOC import opencalphad as oc
#from pyOC import PhaseStatus as phStat
#from pyOC import GridMinimizerStatus as gmStat
from openiec.calculate.coherentenergy import CoherentGibbsEnergy



def test():
    #print(os.getpwd())    
    # Given temperature.
    T= 2773
    #Given pressure
    P=13015
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
    phasenames = ('C1_FCC', 'LIQUID')
    ##set verbolsity
    setverb=False
         
    cge = CoherentGibbsEnergy(T, P, db, comps, x0, phasenames,setverb)
    cge.eqfunc()
   

#    # Molar volumes of pure components to construct corresponding molar volume database.
#    # Molar volume of Al.
#    val = "10.269*10.0**(-6.0) + (3.860*10.0**(-5)*10.0**(-6.0))*T**1.491"
#    # Molar volume of Ni.
#    vni = "6.718*10.0**(-6.0) + (2.936*10.0**(-5)*10.0**(-6.0))*T**1.355"
#    # Molar volume of Cr.
#    vcr = "7.23*10.0**(-6.0)"
#    purevms = [[vni, val, vcr], ]*2
#
#    # A composition range for searching initial interfacial equilirium composition.
#    limit = [0.0001, 0.3]  
#    dx = 0.1
#
#    # Call the module for calculating coherent interfacial energies.
#    sigma = SigmaCoherent(
#        T=T,
#        x0=x0,
#        db=db,
#        comps=comps,
#        phasenames=phasenames,
#        purevms=purevms,
#        limit=limit,
#        dx=dx,
#    )
#
#    # Print the calculated interfacial energy with xarray.Dataset type.
#    print(sigma, "\n")
#    # Print the calculated interfacial energy with xarray.DataArray type.
#    print(sigma.Interfacial_Energy, "\n")
#    # Print the calculated interfacial energy value.
#    print(sigma.Interfacial_Energy.values, "\n")
#
#    # Output
#    """
#    <xarray.Dataset>
#    Dimensions:                     (Components: 3)
#    Coordinates:
#    * Components                  (Components) <U2 'NI' 'AL' 'CR'
#    Data variables:
#        Temperature                 int64 1273
#        Initial_Alloy_Composition   (Components) float64 0.8119 0.18 0.0081
#        Interfacial_Composition     (Components) float64 0.798 0.1931 0.008863
#        Partial_Interfacial_Energy  (Components) float64 0.02378 0.02378 0.02378
#        Interfacial_Energy          float64 0.02378 
#
#    <xarray.DataArray 'Interfacial_Energy' ()>
#    array(0.023783372796539266) 
#
#    0.023783372796539266 
#    """
#
#
if __name__ == "__main__":
    test()
