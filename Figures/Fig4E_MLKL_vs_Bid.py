# Runs ANRM 1.0 (Irvin et. al 2013) under a range of pro-apoptotic and pro-necrotic caspase 8
# concentrations.

import numpy as np
import pylab as p
import matplotlib as mpl
import pickle
import calibratortools as ct
import simulator_1_0 as sim

# ----------Model and Initial Conditions----------------
from anrm.irvin_anrm_model  import model

range_Bid   = [0, 4000, 8000, 12044, 16000, 20000]

#-----------Calibrated Parameters-----------------------
position = pickle.load(open('irvin_anrm_model_fitted_params.pkl'))

#-----------Simulator Settings--------------------------
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,28800,3000) #8hrs converted to seconds (3000 timepoints)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-5
sims.atol = 1e-5

solve = sim.Solver(sims)
solve.run()

delta_td = []
apopt_td = []
necro_td = []

p.ion()
yout = []

condition_variable = 'Bid_0'
graph_name = 'Bid'
observable = 'Obs_MLKL'
graph_obs  = 'MLKL Activation'
rangecv = range_Bid

for i in rangecv:
    #-----------Initial Conditions--------------------------
    ic_params  = model.parameters_initial_conditions()
    conditions = ct.initial_conditions([condition_variable, 'TNFa_0'], [i, 1500], ic_params)
    ysim = solve.simulate(position = position, observables=True, initial_conc = conditions)
    
    #-----------Plot Parp and MLKL--------------------------
    yout.append(ct.extract_records(ysim, [observable]))

for j in range(len(rangecv)):
    p.plot(sims.tspan/3600.0, yout[j], label = '%s %s per cell' % (rangecv[j], graph_name), linewidth = 3)

p.title('%s in WT cells with 25ng/mL TNF and varying initial %s concentrations' % (graph_obs, graph_name))
p.xlabel('time [hrs]')
p.ylabel('%s [molecules per cell]' % graph_obs)
p.legend(bbox_to_anchor = [0.9, -0.25])

