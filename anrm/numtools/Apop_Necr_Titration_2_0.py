# Runs ANRM 1.0 (Irvin et. al 2013) under a range of pro-apoptotic and pro-necrotic caspase 8
# concentrations.

import numpy as np
import pylab as p
import matplotlib as mpl
import pickle
import calibratortools as ct
import simulator_1_0 as sim

# ----------Model and Initial Conditions----------------
from anrm.irvin_mod_v5_tester  import model
#from anrm.irvin_mod_v5_wo_po4bid  import model
#from anrm.irvin_mod_v4_tester  import model

range_proC8 = [10, 4875, 9750, 19500, 39022, 78000] #starting at zero causes NaNs when you normalize the data.
range_C6    = [10, 100, 1000, 10000, 10000]
range_cFlip = [0, 4875, 9750, 19500, 39022, 78000]
range_Bid   = [0, 4000, 8000, 12044, 16000, 20000]
#range_Bid   = [0, 1000, 2000, 8000]
range_RIP1  = [0, 4000, 8000, 12044, 16000, 20000]
range_BidK  = [0, 1000, 2500, 5000, 7500, 10000]
Zinkel_KOs  = [(['TNFa_0'],[1500]),(['Bid_0', 'TNFa_0'],[0, 1500]),(['Bak_0', 'Bax_0', 'TNFa_0'], [0, 0, 1500]), (['Bid_0', 'Bak_0', 'Bax_0', 'TNFa_0'], [0, 0, 0, 1500])]
Zinkel_KO_names = ['WT', '-/-Bid', 'DKO', 'TKO']


#-----------Calibrated Parameters-----------------------
position = pickle.load(open('CompII_Hypthesis_123_newtopology_2run_v40_Position.pkl'))
#position = pickle.load(open('CompII_Hyp_123_Bid_Hyp0_newtopology_1run_v4_Position.pkl'))
#position = pickle.load(open('CompII_Hypthesis_123_addeddata_4run_v41_Position.pkl'))

#-----------Simulator Settings--------------------------
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,43200,3000) #12hrs converted to seconds (3000 timepoints)
#sims.tspan = np.linspace(0,28800,3000) #8hrs converted to seconds (3000 timepoints)
#sims.tspan = np.linspace(0,14400,3000) #4hrs converted to seconds (3000 timepoints)
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

print len(model.parameters_rules())
condition_variable = 'Bid_0'
graph_name = 'Bid'
observable = 'RIP1_FADD'
graph_obs  = 'RIP1 FADD Association'
rangecv = range_Bid

for i in rangecv:
    #-----------Initial Conditions--------------------------
    ic_params  = model.parameters_initial_conditions()
    #conditions = ct.initial_conditions([condition_variable], [i], ic_params)
    #conditions = ct.initial_conditions(i[0], i[1], ic_params) #Use with Zinkel KOs
    conditions = ct.initial_conditions([condition_variable, 'TNFa_0', 'Bak_0', 'Bax_0'], [i, 1500, 0, 0], ic_params)
    ysim = solve.simulate(position = position, observables=True, initial_conc = conditions)
    
    #-----------Plot Parp and MLKL--------------------------
    yout.append(ct.extract_records(ysim, [observable]))

for j in range(len(rangecv)):
    p.plot(sims.tspan/3600.0, yout[j], label = '%s %s per cell' % (rangecv[j], graph_name), linewidth = 3)

yt = yout[0].T
for j in range(len(rangecv))[1:]:
    yt = np.vstack((yt, yout[j].T))

yt = yt.T
for t in range(len(sims.tspan)):
    print '%5.3f, %12.3f, %12.3f, %12.3f, %12.3f, %12.3f, ' % tuple(yt[t])
    #print '%5.3f, %12.3f, %12.3f, %12.3f, ' % tuple(yt[t])
p.title('%s in TKO cells with 25ng/mL TNF and varying initial %s concentrations' % (graph_obs, graph_name))
#p.title('%s in WT, Bid -/-, DKO and TKO cells with 25ng/mL TNF' % (graph_obs))
p.xlabel('time [hrs]')
p.ylabel('%s [molecules per cell]' % graph_obs)
p.legend(bbox_to_anchor = [0.9, -0.25])


"""
#Figure formatting
font_a = {'family':'Helvetica', 'fontweight':'bold', 'size':'12'}
font_b = {'family':'Helvetica', 'fontweight':'normal', 'size':'12'}

#p.title('%s in WT cells with 25ng/mL TNF and varying initial %s concentrations' % (graph_obs, graph_name))
#p.title('%s in WT, Bid -/-, DKO and TKO cells with 25ng/mL TNF' % (graph_obs))
p.xlabel('time [hrs]', **font_b)
p.ylabel('%s [molecules per cell]' % graph_obs, **font_b)
p.xticks(**font_a)
p.yticks(**font_a)
p.legend(bbox_to_anchor = [1.1, -0.25])
[i.set_linewidth(3) for i in ax.spines.itervalues()]"""