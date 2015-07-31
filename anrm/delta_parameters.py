# Fits ANRM 1.0 (Irvin et. al 2013) against single-cell measurements
# of caspase reporters.

import pickle
import bayessb
import numpy as np
from anrm.numtools import simulator_1_0 as sim

from anrm.irvin_anrm_wo_bidpo4 import model

#-----------Previously Calibrated Parameters------------
position = pickle.load(open('irvin_anrm_wo_bidpo4_fitted_params.pkl'))

#----Simulator Settings----
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,36000,1000) #10hrs converted to seconds (1000 timepoints)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-5
sims.atol = 1e-5

solve = sim.Solver(sims)
solve.run()

# print some information about the maximum-likelihood estimate parameter set
print
print '%-10s %-12s %-12s %s' % ('parameter', 'actual', 'fitted', 'log10(fit/actual)')
fitted_values = solve.cur_params(position)[solve.estimate_idx]
for param, new_value in zip(sims.estimate_params, fitted_values):
    change = np.log10(new_value / param.value)
    values = (param.name, param.value, new_value, change)
    print '%-10s, %-12.2g, %-12.2g, %-+6.2f,' % values

