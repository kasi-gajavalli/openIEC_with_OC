"""
Obtain quantities correlating with thremodynamic equilibrium calculation using the pycalphad package.
"""

import matplotlib.pyplot as plt
from pycalphad import equilibrium
from pycalphad import Database, Model
import pycalphad.variables as v
import numpy as np
import math

"""
Obtain quantities correlating with thremodynamic equilibrium calculation using the pyOC package.
"""
import sys
sys.path.append('/home/kg245220/code/pyOC')
import pyOC
import numpy as np
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
from pyOC import GridMinimizerStatus as gmstat


class CoherentGibbsEnergy_OC(object): 
    def __init__(self, T, P, db, comps, x0, phasename,setverb,pathname):
        self.T = T
        self.P = P
        self.db = db
        self.comps = comps
        self.phasename = phasename
        
    def eqfunc(self, x):
        """Calculate Phase Equilibrium"""
        # setting verbosity (True or False - default), if set to yes, in particular, when getters are called the returned values are displayed in a comprehensive way
        oc.setVerbosity(self.setverb)
        
         # tdb filepath
        tdbFile=self.pathname+self.db
        #print(tdbFile)
         # reading tdb
        oc.readtdb(tdbFile,self.comps)
        oc.setPhasesStatus((self.phasename[0],self.phasename[1]),phStat.Entered)
       
        #Equilibrium 
        variable = oc.setElementMolarAmounts(self.x)+ [oc.setTemperature(self.T), oc.setPressure(self.P)]
        xs = [v.X(each) for each in self.comps if each != "VA"]
        xxs = [xs[i] for i in range(1, len(xs))]
        xxxs = xxs + [v.T, v.P]
        var = {xxxs[i]: variable[i] for i in range(len(variable))}
        eq_result = oc.calculateEquilibrium(self.db, self.comps, self.phasename, var)
        return eq_result
    
    def phase(self, x, **kwargs):
        """
        The string name of the phase in equilibrium at the conditions.
        """
        phasearray = self.eqfunc(x).Phase.sel(P=self.P, T=self.T, **kwargs)
        return phasearray.values.flatten()
        
        phasesAtEquilibrium=oc.getPhasesAtEquilibrium()
        
        
        self.pec = phasesAtEquilibrium.getPhaseElementComposition()
        self.pcc = phasesAtEquilibrium.getPhaseConstituentComposition() 
        self.ps = phasesAtEquilibrium.getPhaseSites()
        self.ma = phasesAtEquilibrium.getPhaseMolarAmounts()        
        self.mu = oc.getChemicalPotentials()
        self.cd = oc.getConstituentsDescription()
        self.mass=oc.getScalarResult('B')
        
    def getGibbsEnergy(self):
        return self.gibbs
         
    def getMolarAmounts(self):
        return self.ma
    
    def getPhaseElementComposition(self):
        return self.pec
    
    def getPhaseConstituentComposition(self):
        ###Mole fractions of components for phases in equilibrium.
        return self.pcc
       
    def getPhaseSites(self):
        return self.ps
        
        
    def getChemicalPotential(self):
        #Chemical potentials of components in bulk phases
        return self.mu
        
    def getConstituentsDescription(self):
        return self.cd
    
    def getMass(self):
         return self.mass
     
    def chemicalpotential(self, x, **kwargs):
        componentname = [each for each in self.comps if each != "VA"]
        chemicalpotential = (
            self.eqfunc(x).MU.sel(P=self.P, T=self.T, component=componentname, **kwargs)
        ).values
        return chemicalpotential.flatten()