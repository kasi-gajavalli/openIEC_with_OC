

"""
Obtain quantities correlating with thremodynamic equilibrium calculation using the pyOC package.
"""
#import sys
#sys.path.append('/home/kg245220/code/pyOC')
import pyOC
#import numpy as np
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
from pyOC import GridMinimizerStatus as gmstat


class CoherentGibbsEnergy_PyOC(object):
    
    def __init__(self, T, P, db, comps, phasenames, setverb):
        self.T = T
        self.P = P
        self.db = db
        self.comps = comps
        self.phasenames = phasenames
        self.setverb = setverb
        self.eq_val = None
        
    def readDatabase(self):
        oc.readtdb(self.db,self.comps)
        
    def suspendPhases(self, suspended_phase):
        
        oc.setPhasesStatus((suspended_phase,),phStat.Suspended)
        
        entered_phase = self.phasenames[self.phasenames==suspended_phase]
        
        oc.setPhasesStatus((entered_phase,),phStat.Entered)
    
    def multiPhaseAnalysis(self):
        oc.setPhasesStatus((self.phasenames[0],self.phasenames[1]),phStat.Entered)
        
    def eqfunc(self, x, calc_value):
        
        # setting verbosity (True or False - default), if set to yes, in particular, when getters are called the returned values are displayed in a comprehensive way
        oc.setVerbosity(self.setverb)
        
         # set pressure
        oc.setPressure(self.P)
        
          # set temperature
        oc.setTemperature(self.T)
    
        # set initial molar amounts
        oc.setElementMolarAmounts(x)
       
        #Equilibrium 
        oc.calculateEquilibrium(gmstat.Off)
        
        if calc_value == 'gibbs':
            self.eq_val = oc.getGibbsEnergy()
        
        elif calc_value == 'mu':
            self.eq_val = oc.getChemicalPotentials()
        
        elif calc_value == 'cd':
            self.eq_val = oc.getConstituentsDescription()
        
        elif calc_value == 'mass':
            self.eq_val = oc.getScalarResult('B')
        
        elif calc_value == 'pec':
            self.eq_val = oc.getPhasesAtEquilibrium().getPhaseElementComposition()
        
        elif calc_value == 'pcc':
            self.eq_val = oc.getPhasesAtEquilibrium().getPhaseConstituentComposition() 
            
        elif calc_value == 'ps':
            self.eq_val = oc.getPhasesAtEquilibrium().getPhaseSites()
        
        elif calc_value == 'ma':
            self.eq_val = oc.getPhasesAtEquilibrium().getPhaseMolarAmounts()
        
        return self.eq_val
        
    def Gibbsenergy(self, x):
        gibbs = self.eqfunc(x, 'gibbs')
        return gibbs
    
    def MolarAmounts_PyOC(self, x):
        ma = self.eqfunc(x, 'ma')
        return ma
    
    def PhaseElementComposition_PyOC(self, x):
        pec = self.eqfunc(x, 'pec')
        return pec
    
    def PhaseConstituentComposition_PyOC(self, x):
        pcc = self.eqfunc(x, 'pcc')
        return pcc
       
    def phasesites(self, x):
        ps = self.eqfunc(x, 'ps')
        return ps
        
    def chemicalpotential(self, x):
        mu = self.eqfunc(x, 'mu') # its a dictionary converted to the list of values
        return list(mu.values())
        
    def ConstituentsDescription_PyOC(self, x):
        cd = self.eqfunc(x, 'cd')
        return cd
    
    def mass(self, x):
        mass = self.eqfunc(x, 'mass')
        return mass
     
#   def eq_func_with_missicibilityGap(self, gmStat, x, calc_value):
#        
#        # setting verbosity (True or False - default), if set to yes, in particular, when getters are called the returned values are displayed in a comprehensive way
#        oc.setVerbosity(self.setverb)
#        
#         # set pressure
#        oc.setPressure(self.P)
#        
#          # set temperature
#        oc.setTemperature(self.T)
#    
#        # set initial molar amounts
#        oc.setElementMolarAmounts(x)
#       
#        #Equilibrium 
#        if gmStat=='Off':
#            oc.calculateEquilibrium(gmstat.Off)
#        elif gmStat=='On':
#            oc.calculateEquilibrium(gmstat.On)
#        else:
#            raise ValueError('No suitable parameter for gmstat found: Choose from On/Off')
#        
#        #Equilibrium 
#        oc.calculateEquilibrium(gmstat.Off) ##### change afterwards
#        
#        if calc_value == 'gibbs':
#            self.eq_val = oc.getGibbsEnergy()
#        
#        elif calc_value == 'mu':
#            self.eq_val = oc.getChemicalPotentials()
#        
#        elif calc_value == 'cd':
#            self.eq_val = oc.getConstituentsDescription()
#        
#        elif calc_value == 'mass':
#            self.eq_val = oc.getScalarResult('B')
#        
#        elif calc_value == 'pec':
#            self.eq_val = oc.getPhasesAtEquilibrium().getPhaseElementComposition()
#        
#        elif calc_value == 'pcc':
#            self.eq_val = oc.getPhasesAtEquilibrium().getPhaseConstituentComposition() 
#            
#        elif calc_value == 'ps':
#            self.eq_val = oc.getPhasesAtEquilibrium().getPhaseSites()
#        
#        elif calc_value == 'ma':
#            self.eq_val = oc.getPhasesAtEquilibrium().getPhaseMolarAmounts()
#        
#        return self.eq_val
        


