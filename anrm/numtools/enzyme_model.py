from pysb import *
from pysb.util import alias_model_components

import Enzyme_Kinetics_Example as ee

Model('m')

ee.enzyme_monomers()
ee.enzyme_initials()
ee.enzyme_kinetics()
#ee.observables()

alias_model_components()