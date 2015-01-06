# Fits ANRM 1.0 (Irvin et. al 2013) against single-cell measurements
# of caspase reporters.

import pickle
import scipy
import numpy as np
import simulator_1_0 as sim
import matplotlib.pyplot as plt
from SloppyCell.ReactionNetworks import *

###from anrm.irvin_mod_v5_tester import model
from anrm.irvin_mod_v5_wo_po4bid import model

#-----------Previously Calibrated Parameters------------
###position = pickle.load(open('CompII_Hypthesis_123_newtopology_1run_v4_Position.pkl'))
position = pickle.load(open('CompII_Hyp_123_Bid_Hyp0_newtopology_1run_v4_Position.pkl'))

#-----------Import Hessian------------------------------
###h = pickle.load( open('pysb_hessian_CompII_Hyp_123_Apop2_Necr1.pkl'))
h = pickle.load( open('pysb_hessian_CompII_Hyp_123_Bid_Hyp_0_Apop2_Necr1.pkl'))

###h = pickle.load( open('sc_hessian_necro'))

#-----------Plot Eigen-spectra--------------------------
e, v = Utility.eig(h)
e = scipy.real(e)
    
width = (234/len(e))**0.25 * 0.25
l = Plotting.plot_eigval_spectrum(e, offset = 0.15,
                                      widths=0.7, lw=width)

#-----------Print Largest Eigenvalues and Eigenvectors--
"""retrieving parameters and values"""
#----Simulator Settings----
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,36000,1000) #10hrs converted to seconds (1000 timepoints)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-5
sims.atol = 1e-5

solve = sim.Solver(sims)
solve.run()

fig = Plotting.figure()

for i in [0, 1, 2, 3, 4, 5, 6]:
    print "the %sthlargest eigenvalue =%s " % (i+1, e[i])
    print
    print '%-30s %-12s %s' % ('parameter', 'value', 'eigenvector_component')

    fitted_values = solve.cur_params(position)[solve.estimate_idx]
    for param, value, ev_component in zip(sims.estimate_params, fitted_values, v[i]):
        values = (param.name, value, ev_component)
        print '%-30s %-12.2g %-+6.2g' % values

Plotting.show()
