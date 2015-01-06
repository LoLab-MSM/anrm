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
from pysb.macros import *

Model()
# SECTION ONE: Fas Signalling
# ===========================
# This contains FasL binding and assembly of the DISC

def enzyme_monomers():
    
    """ Basic single substrate enzyme kinetics.
        A + B <> A:B >> C
    """
    Monomer('A', ['a', 'b'])
    Monomer('B', ['c', 'd'])
    Monomer('C', ['e', 'f'])
    Monomer('D', ['g', 'h'])
    alias_model_components()

def enzyme_initials():
    Parameter('A_0' ,  60000) # 6000 corresponds to 100ng/ml TNFa
    Parameter('B_0' ,  8000) # 4800 receptors per cell
    Parameter('C_0' ,  60000) # 6000 corresponds to 100ng/ml TNFa
    Parameter('D_0' ,  8000) # 4800 receptors per cell
    alias_model_components()
    
    Initial(A(a=None, b=None), A_0)
    Initial(B(c=None, d = None), B_0)
    Initial(C(e=None, f=None), C_0)
    Initial(D(g=None, h = None), D_0)

"""def enzyme_kinetics():
    Rule('rxn1', A(bf=None) + A(bf=None) <> A(bf = 1)%A(bf = 1), Parameter('k1', 0.1), Parameter('k2', 1e-3))
    Rule('rxn2', A(bf=1)%A(bf=1) >> A(bf=None) + C(bf=None), Parameter('k3',1e-3))
    alias_model_components ()"""
def enzyme_kinetics():
    #bind_complex(A(a=ANY, b=None), 'b', C(e=2)%D(g=2), 'h', [1e-4, 1e-1])
    catalyze_complex(C(e=2)%D(g=2), 'h', A(a=1, b=None)%B(c=1, d=None), 'b', A(a=1, b=None)%B(c=1, d=2)%B(c=2, d=None), [1e-4, 1e-1, 1e-1])
    
    alias_model_components()
         
def observables():
    Observable('obs_A', A(a = None))
    Observable('obs_AB', B(c = None))
    Observable('obs_C', C(e = None))
    alias_model_components()

