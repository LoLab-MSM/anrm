
from pylab import *
from scipy import *
from SloppyCell.ReactionNetworks import *

from anrm.irvin_mod_v5_tester import model

import pickle
import numpy as np
import scipy as sp
import simulator_1_0 as sim
import sys

sys.setrecursionlimit(5000)
experiment_list = ['Necr1']

# -------SETTING UP NETWORK-------
# import network topology
net = IO.from_SBML_file('irvin_mod_v5_tester.sbml', experiment_list[0])

import ANRM_expt

# set initial conditions
for i in ANRM_expt.init_conditions[experiment_list[0]].keys():
    net.set_var_ic(i, ANRM_expt.init_conditions[experiment_list[0]][i])

# set observables
for i in experiment_list:
    net.add_species(ANRM_expt.obs[i][0], 0)
    net.add_assignment_rule(ANRM_expt.obs[i][0], ANRM_expt.obs[i][1])

"""Set Parameters"""

# import position and assign parameter values
position = pickle.load(open('CompII_Hypthesis_123_newtopology_1run_v4_Position.pkl'))

# use methods in simulator_1_0 to convert parameters to a set that is
# compatible with sloppy cell
sims = sim.Settings()
sims.model = model
sims.estimate_params = model.parameters_rules()

solve = sim.Solver(sims)
solve.estimate_idx =[i for i, p in enumerate(model.parameters) if p in sims.estimate_params]

position_values = solve.cur_params(position)[solve.estimate_idx]

params_list = zip(sims.estimate_params.keys(), position_values)

params = KeyedList(params_list)

# time course data
m = Model([ANRM_expt.Necr1],[net])

print 'Initial cost:', m.cost(params)

for id, val in params.items():
    m.AddResidual(Residuals.PriorInLog('prior_on_%s' % id, id, sp.log(val),
                                       sp.log(10)))

