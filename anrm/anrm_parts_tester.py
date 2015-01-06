from pysb import *
from pysb.util import alias_model_components
from earm import shared

from irvin_modules_v6 import *

Model('m')

TNFa_to_ComplexI_Monomers()
TNFa_to_ComplexI_Initials()
TNFa_to_ComplexI()