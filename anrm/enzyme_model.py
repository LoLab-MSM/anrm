from pysb import *
from pysb.util import alias_model_components
from earm import shared

from Enzyme_Kinetics_Example import *

Model('m')

enzyme_monomers()
enzyme_initials()
enzyme_kinetics()
observables()