

"""
Obtain quantities correlating with thremodynamic equilibrium calculation using the pyOC package.
"""
import pyOC
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
from pyOC import GridMinimizerStatus as gmstat


class CoherentGibbsEnergy_OC(object):
    
    def __init__(self):
        self.pec = None
        self.pcc = None
        self.ps = None
        self.ma = None
        self.cd = None
        self.gibbs = 0.
        self.mu = 0.
        self.mass = 0.
        
    def readDatabase(self, pathname, db, comps):
        oc.readtdb(pathname+db,comps)
        
    def suspendPhases(self, suspended_phase, phasenames):
        
        oc.setPhasesStatus((suspended_phase,),phStat.Suspended)
        
        entered_phase = phasenames[phasenames==suspended_phase]
        
        oc.setPhasesStatus((entered_phase,),phStat.Entered)
    
    def multiPhaseAnalysis(self, phasenames):
        oc.setPhasesStatus((phasenames[0],phasenames[1]),phStat.Entered)
        
    def calculateEquilibria(self, gmStat, setverb, T, P, x0):
        
        # setting verbosity (True or False - default), if set to yes, in particular, when getters are called the returned values are displayed in a comprehensive way
        oc.setVerbosity(setverb)
        
         # set pressure
        oc.setPressure(P)
        
          # set temperature
        oc.setTemperature(T)
    
        # set initial molar amounts
        oc.setElementMolarAmounts(x0)
       
        #Equilibrium 
        if gmStat=='Off':
            oc.calculateEquilibrium(gmstat.Off)
        elif gmStat=='On':
            oc.calculateEquilibrium(gmstat.On)
        else:
            raise ValueError('No suitable parameter for gmstat found: Choose from On/Off')
        
        self.gibbs=oc.getGibbsEnergy()
        self.mu = oc.getChemicalPotentials()
        self.cd = oc.getConstituentsDescription()
        self.mass=oc.getScalarResult('B')
        
        phasesAtEquilibrium=oc.getPhasesAtEquilibrium()
        self.pec = phasesAtEquilibrium.getPhaseElementComposition()
        self.pcc = phasesAtEquilibrium.getPhaseConstituentComposition() 
        self.ps = phasesAtEquilibrium.getPhaseSites()
        self.ma = phasesAtEquilibrium.getPhaseMolarAmounts()
        
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
     
   


