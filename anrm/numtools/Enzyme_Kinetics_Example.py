"""
    Overview
    ========
    
    Provides a model of enzyme kinetics. And a second reaction model which does not
    describe enzyme kinetics.
    --
    """

import numpy
from pysb import *
from pysb.util import alias_model_components
from pysb.bng import *
from pysb.macros import *
from pysb.integrate import odesolve
# SECTION ONE: Fas Signalling
# ===========================
# This contains FasL binding and assembly of the DISC

def enzyme_monomers():
    """
    Basic single substrate enzyme kinetics.
        A + B <> A:B >> C"""
    Monomer('A'  ,  ['bf'])
    Monomer('B' ,   ['bf'])
    Monomer('C' ,   ['bf'])
    alias_model_components()

def enzyme_initials():
    Parameter('A_0' ,  600000) # 6000 corresponds to 100ng/ml TNFa
    Parameter('B_0' ,  8000) # 4800 receptors per cell
    alias_model_components()
    
    Initial(A(bf=None), A_0)
    Initial(B(bf=None), B_0)

"""def enzyme_kinetics():
    Rule('rxn1', A(bf=None) + A(bf=None) <> A(bf = 1)%A(bf = 1), Parameter('k1', 0.1), Parameter('k2', 1e-3))
    Rule('rxn2', A(bf=1)%A(bf=1) >> A(bf=None) + C(bf=None), Parameter('k3',1e-3))
    alias_model_components ()"""
def enzyme_kinetics():
    bind(A(bf=None), 'bf', B(bf=None), 'bf', [0.1, 1e-3])
    Rule('rxn2', A(bf=1)%A(bf=1) >> A(bf=None) + C(bf=None), Parameter('k3',1e-3))
    alias_model_components()

def observables():
    Observable('obs_A', A(bf = None))
    Observable('obs_AB', A(bf = 1)%B(bf=1))
    Observable('obs_C', C(bf = None))
    alias_model_components()


"""

# Simulate the model through 40 seconds
time = numpy.linspace(0, 40, 100)
print "Simulating..."
x = odesolve(model, time)"""



