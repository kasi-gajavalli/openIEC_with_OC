#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 13:02:27 2020

@author: kg245220
"""

"""According to given conditions and input parameters, calculate corresponding interfacial energies.
"""
#import sys
#sys.path.append('/home/kg245220/code/openiec_with_OC/openiec/property')

from pycalphad import Database
from openiec.property.PyOC_coherentenergy import CoherentGibbsEnergy_PyOC
#sys.path.append('/home/kg245220/code/openiec_with_OC/openiec/model')
from openiec.model.sigmacoint import SigmaCoherentInterface
from openiec.property.molarvolume import MolarVolume, InterficialMolarVolume
#sys.path.append('/home/kg245220/code/openiec_with_OC/openiec/calculate')
from openiec.calculate.minimize import SearchEquilibrium, ComputeEquilibrium
from openiec.utils.decorafunc import wraptem
import numpy as np
from xarray import Dataset



def SigmaCoherent(T, x0, db, comps, phasenames, purevms, intervms=[], limit=[0, 1.0], dx=0.01, P=1E5, setverb=False):
    """
    Calculate the coherent interfacial energy in alloys.

    Parameters
    -----------
    T: float
        Given temperatur
    x0: list
        Initial alloy composition.
    db : Database
        Database containing the relevant parameters.
    comps : list
        Names of components to consider in the calculation.
    phasenames : list
        Names of phase model to build.    
    limit: list
        The limit of composition for searching interfacial composition in equilibrium.
    purevms: list
        The molar volumes of pure components.
    dx: float
        The step of composition for searching interfacial composition in equilibrium.

    Returns:   
    -----------
    Components：list of str
        Given components.
    Temperature: float
        Given temperature.
    Initial_Alloy_Composition: list
        Given initial alloy composition.
    Interfacial_Composition: list
        Interfacial composition of the grid minimization.
    Partial_Interfacial_Energies: list
        Partial interfacial energies of components.
    Interfacial_Energy: float    
        Requested interfacial energies.

    Return type: xarray Dataset
    """
    phasevm = [MolarVolume(Database(db), phasenames[i], list(comps), purevms[i]) for i in range(2)] #### ASK ???????????????????????
    
    _vmis = InterficialMolarVolume(*phasevm)

    """decorate the _vmis to release the constains on temperature"""
    vmis = [wraptem(T, f) for f in _vmis]

    
    model = CoherentGibbsEnergy_PyOC(T, P, db, comps, phasenames, setverb)
    model.readDatabase()
    model.multiPhaseAnalysis()
    
    """Chemical potentials in two-phase equilibrium"""
    mueq = model.chemicalpotential(x0)

    """Chemical potentials in two bulk phases"""
    model_phase = [
        CoherentGibbsEnergy_PyOC(T, P, db, comps, phasenames[i], setverb) for i in range(len(phasenames))
    ]
    alphafuncs, betafuncs = [each.chemicalpotential for each in model_phase]

    sigma_model = SigmaCoherentInterface(alphafuncs, betafuncs, mueq, vmis, comps)

    components = [each for each in comps if each != "VA"]
    cum = int(len(components) - 1)
    print(
        "\n******************************************************************************\nOpenIEC is looking for interfacial equilibirium coposition.\nFor more information visit http://....../openiec."
    )
    x_s = SearchEquilibrium(sigma_model.objective, [limit] * cum, [dx] * cum)
    x_c = ComputeEquilibrium(sigma_model.objective, x_s["x"])
    print(
        "******************************************************************************\n\n"
    )
    sigma = sigma_model.infenergy(x_c)
    print(sigma)

    xx0 = [1.0 - sum(list(x0.values()))] + list(x0.values())
    xx_c = [1.0 - sum(list(x_c))] + list(x_c)
    sigmapartial = list(np.array(sigma).flatten())
    sigmaavg = np.average(sigma)

    res = Dataset(
        {
            "Components": components,
            "Temperature": T,
            "Initial_Alloy_Composition": ("Components", xx0[1:]),
            "Interfacial_Composition": ("Components", xx_c),
            "Partial_Interfacial_Energy": ("Components", sigmapartial),
            "Interfacial_Energy": sigmaavg
        }
    )

    return res

