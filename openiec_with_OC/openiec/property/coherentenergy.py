"""
Obtain quantities correlating with thremodynamic equilibrium calculation using the pycalphad package.
"""
import sys
sys.path.append('/home/kg245220/code/pyOC')
import pyOC
import numpy as np
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
from pyOC import GridMinimizerStatus as gmstat


class CoherentGibbsEnergy(object): 
    def __init__(self, T, P, db, comps, x0, phasename,setverb,pathname):
        self.T = T
        self.P = P
        self.db = db
        self.x0 = x0
        self.comps = comps
        self.phasename = phasename
        self.setverb = setverb
        self.pathname=pathname
        self.pef = None
        self.pcc = None
        self.ps = None
        self.ma = None
        self.gibbs = 0.
        self.mu = 0.
        self.cd = None
        self.mass=0.
        
    
    def eqfunc(self,gmStat,element,multi_phase_analysis):
        
        # setting verbosity (True or False - default), if set to yes, in particular, when getters are called the returned values are displayed in a comprehensive way
        oc.setVerbosity(self.setverb)
        
         # tdb filepath
        tdbFile=self.pathname+self.db
        #print(tdbFile)
         # reading tdb
        oc.readtdb(tdbFile,self.comps)
        
        if multi_phase_analysis == False:
            if element == 'O':
                oc.setPhasesStatus((self.phasename[1],),phStat.Suspended)
                oc.setPhasesStatus((self.phasename[0],),phStat.Entered)
            else:
                oc.setPhasesStatus(('*',),phStat.Suspended)
                oc.setPhasesStatus((self.phasename[0],),phStat.Entered)
                
        else:
            oc.setPhasesStatus((self.phasename[0],self.phasename[1]),phStat.Entered)
        
         # set pressure
        oc.setPressure(self.P)
        
        
          # set temperature
        oc.setTemperature(self.T)
    
        # x0=x0
        oc.setElementMolarAmounts(self.x0)
       
        #Equilibrium 
        
        if gmStat=='Off':
            oc.calculateEquilibrium(gmstat.Off)
        elif gmStat=='On':
            oc.calculateEquilibrium(gmstat.On)
        else:
            raise ValueError('No suitable parameter for gmstat found: Choose from On/Off')
        
        self.gibbs=oc.getGibbsEnergy()
        
        
        phasesAtEquilibrium=oc.getPhasesAtEquilibrium()
        
        
        self.pef = phasesAtEquilibrium.getPhaseElementComposition()
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
    
    def getPhaseElementFraction(self):
        return self.pef
    
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

