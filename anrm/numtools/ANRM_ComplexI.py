from pysb import *
from pysb.util import alias_model_components

import ComplexI_partial_rules as irvin

Model('model')

irvin.TNFa_to_ComplexI_Monomers()

irvin.TNFa_to_ComplexI_Initials()
irvin.TNFa_to_ComplexI()

alias_model_components()