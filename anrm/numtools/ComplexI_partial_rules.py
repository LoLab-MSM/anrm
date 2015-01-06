"""
    Overview
    ========
    
    Provides a model of enzyme kinetics. And a second reaction model which does not
    describe enzyme kinetics.
    --
    """
import numpy
from pysb.macros import *

import numpy
from pysb import *
from pysb.util import alias_model_components
from pysb.bng import *
from pysb.integrate import odesolve
# SECTION ONE: Fas Signalling
# ===========================
# This contains FasL binding and assembly of the DISC

def TNFa_to_ComplexI_Monomers():
    
    """ Declares TNFa, TNFR1, TRADD, RIP1, TRAF2, IAP, NKS (NFkB signaling complex),
        NFkB, CYLD and FADD. Binding sites are named after (when known) subdomains that
        facilitate binding. Otherwise the binding site will reflect what it binds with.
        
        The 'state' site denotes various localization and/or activity states of a
        Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
        localization.
        """
    Monomer('TNFa'  ,   ['brec'])           #TNFa
    Monomer('TNFR1' ,   ['blig', 'bDD'])    #TNFR1
    Monomer('TRADD' ,   ['bDD1', 'bDD2'])   #TRADD
    Monomer('RIP1'  ,   ['bDD', 'btraf', 'bRHIM', 'bMLKL', 'state'], {'state':['unmod', 'ub', 'po4', 'trunc']})   #RIP1
    Monomer('TRAF'  ,   ['brip', 'bciap', 'zfinger'])   #TRAF2
    Monomer('cIAP'  ,   ['btraf'])          #cIAP 1 and 2
    Monomer('NSC'   ,   ['bnfkb'])          #Recruitement of LUBAC commits complex I to the NFkB activation pathway. This pathway is approximated by a single activation step which is carried out by a NFkB signaling complex (NSC).
    Monomer('NFkB'  ,   ['bnsc','state'], {'state':['I', 'A']})       #NFkB
    Monomer('CYLD'  ,   ['btraf', 'state'], {'state':['U','T']})  #CYLD
    Monomer('FADD', ['bDD', 'bDED1','bDED2'])    #FADD
    alias_model_components()

def TNFa_to_ComplexI_Initials():
    Parameter('TNFa_0'  ,  600) # 6000 corresponds to 100ng/ml TNFa
    Parameter('TNFR1_0' ,  4800) # 4800 receptors per cell
    Parameter('TRADD_0' ,  9000) # molecules per cell (arbitrarily assigned) 9000
    Parameter('RIP1_0'  , 12044) # molecules per cell (arbitrarily assigned) 12044
    Parameter('TRAF_0'  ,  9000) # molecules per cell (arbitrarily assigned) 9000
    Parameter('cIAP_0'  ,  9000) # molecules per cell (arbitrarily assigned) 9000
    Parameter('NSC_0'   ,     0) # complexes per cell
    Parameter('NFkB_0'    , 50000) # molecules per cell (arbitrarily assigned) 50000
    Parameter('CYLD_0'  ,  9000) # molecules per cell
    Parameter('FADD_0'  ,   8030) # 8300 molecules per cell (J Immunol 2005)
    alias_model_components()
    
    Initial(TNFa(brec=None), TNFa_0)
    Initial(TNFR1(blig=None, bDD=None), TNFR1_0)
    Initial(TRADD(bDD1=None, bDD2=None), TRADD_0)
    Initial(RIP1(bDD=None, btraf=None, bRHIM = None, bMLKL = None, state = 'unmod'), RIP1_0)
    Initial(TRAF(brip=None, bciap=None, zfinger=None), TRAF_0)
    Initial(cIAP(btraf=None), cIAP_0)
    Initial(NSC(bnfkb=None), NSC_0)
    Initial(NFkB(bnsc=None, state = 'I'), NFkB_0)
    Initial(CYLD(btraf=None, state = 'U'),CYLD_0)
    Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0)

def TNFa_to_ComplexI():
    """Reaction network that produces Complex I and activates NFkB [1]. Recruitment of LUBAC and NFkB pathway components to complex I is approximated by a one-step conversion to "NFkB Signaling Complex" (reaction 7). A20 ubiquitylates RIP1 and targets it for degradation. The action of A20 is apporixmated at a single reaction (reaction 8).
        
        1. TNFa + TNFR1 <> TNFa:TNFR1
        2. TNFa:TNFR1 + TRADD <> TNFa:TNFR1:TRADD
        3. TNFa:TNFR1:TRADD + RIP1(unmod) <> TNFa:TNFR1:TRADD:RIP1(unmod)
        4. TNFa:TNFR1:TRADD:RIP1(unmod) + TRAF2 <> TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2
        5. TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2 + cIAP <> TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2:cIAP
        
        6. TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2:cIAP >>TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2:cIAP
        7. TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2 >> NFkB Signaling Complex
        8. TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2 >> TNFa:TNFR1:TRADD + TRAF2
        9. TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2 + CYLD <> TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2:CYLD
        10.TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2:CYLD >> TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2:CYLD
        
        1. Olivier Micheau, Jurg Tschopp, Induction of TNF Receptor I-Mediated Apoptosis via Two Sequential Signaling Complexes, Cell, Volume 114, Issue 2, 25 July 2003, Pages 181-190
        """
    bind(TNFa(brec = None),'brec', TNFR1(blig = None), 'blig', [1e-6, 1e-3])
    bind(TNFR1(blig = ANY, bDD = None), 'bDD', TRADD(bDD1 = None, bDD2 = None), 'bDD1', [1e-6, 1e-3])
    #we need to get bind_complex working!
    Rule('RIP1_to_complex1', TNFR1(blig = ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = None) + RIP1(bDD = None, btraf = None, state =  'unmod') <> TNFR1(blig = ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = 1)%RIP1(bDD = 1, btraf = None, state =  'unmod'), Parameter('k1', 1e-6),Parameter('k2',1e-3))
    Rule('TRAF_to_complex1', TNFR1(blig =  ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)%RIP1(bDD = ANY, btraf = None, state =  'unmod') + TRAF(brip=None) <> TNFR1(blig =  ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)%RIP1(bDD = ANY, btraf = 1, state =  'unmod')%TRAF(brip=1), Parameter('k3', 1e-6), Parameter('k4',1e-3))
    bind(TRAF(bciap =  None), 'bciap', cIAP(btraf = None), 'btraf', [1e-6, 1e-3])
    bind(TRAF(zfinger =  None), 'zfinger', CYLD(btraf = None), 'btraf', [1e-6, 1e-3])

    ComplexI = TNFa(brec = ANY)%TNFR1(blig =  ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)

    Rule('RIP1_Ubiquitination', ComplexI %RIP1(bDD = ANY, btraf = ANY, state =  'unmod')%TRAF(brip = ANY, bciap = ANY) >> ComplexI %RIP1(bDD = ANY, btraf = None, state =  'ub')%TRAF(brip = ANY, bciap = ANY), Parameter('RIP1_ubiq_kc', 1e-1))
    Rule('RIP1_Deubiquitination', ComplexI%RIP1(bDD = ANY, btraf = ANY, state =  'ub')%TRAF(brip = ANY, zfinger = ANY) >> ComplexI%RIP1(bDD = ANY, btraf = None, state =  'unmod')%TRAF(brip = ANY, zfinger = ANY), Parameter('RIP1_deub_kc', 1e-1))
    #Rule('Establish_NFkB_Signaling_Complex', ComplexI%RIP1(bDD = ANY, btraf = ANY, state =  'ub')%TRAF(brip = ANY, bciap = ANY) >> NSC(bnfkb=None), Parameter('NSC_esta_kc', 1e-8))
    #Rule('RIP1_Degradation', ComplexI%RIP1(bDD = ANY, btraf = ANY, state =  'ub')%TRAF(brip = ANY) >> ComplexI + TRAF(brip = None), Parameter('RIP1_degr_kc', 1e-1))




