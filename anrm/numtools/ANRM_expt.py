import numpy as np

from SloppyCell.ReactionNetworks import *

#----Generating Timecourse Data-----
def ydata_fn():
    """return and array of synthesized experimental data. The array is loosely based on published experiments"""
    Apop1_td = 6.0 #six hours
    Apop2_td = 4.0 #four hours
    Necr1_td = 4.0 #four hours

    switchtime_CytoC = 1.0 # [hrs]
    switchtime_cPARP = 0.5 #one hour
    switchtime_MLKL = 1.0 # [hrs]

    Apop1_obs = ['Obs_CytoC'] #Zhang et al. Monitored CytoC (Obs_CytoC) but CytoC does not have switch behavior.
    Apop2_obs = ['Obs_cPARP']
    Necr1_obs = ['Obs_MLKL']
    
    ydata = {}
    ydata['Apop1'] = [np.array([[(Apop1_td-2*switchtime_CytoC),(Apop1_td-switchtime_CytoC), (Apop1_td-switchtime_CytoC/2), (Apop1_td-switchtime_CytoC/4), (Apop1_td-switchtime_CytoC/8), Apop1_td, (Apop1_td+switchtime_CytoC/8), (Apop1_td+switchtime_CytoC/4), (Apop1_td+switchtime_CytoC/2), (Apop1_td+switchtime_CytoC)], [0, 0, 0.05, 0.205, 0.340, 0.5, 0.659, 0.794, 0.95, 1],[0.025, 0.025, 0.025, 0.05, 0.075, 0.1, 0.085, 0.065, 0.05, 0.025]]).T, Apop1_obs]
    ydata['Apop2'] = [np.array([[(Apop2_td-2*switchtime_cPARP),(Apop2_td-switchtime_cPARP), (Apop2_td-switchtime_cPARP/2), (Apop2_td-switchtime_cPARP/4), (Apop2_td-switchtime_cPARP/8), Apop2_td, (Apop2_td+switchtime_cPARP/8), (Apop2_td+switchtime_cPARP/4), (Apop2_td+switchtime_cPARP/2), (Apop2_td+switchtime_cPARP)], [0, 0, 0.05, 0.205, 0.340, 0.5, 0.659, 0.794, 0.95, 1], [0.025, 0.025, 0.025, 0.05, 0.075, 0.1, 0.085, 0.065, 0.05, 0.025]]).T, Apop2_obs]
    ydata['Necr1'] = [np.array([[(Necr1_td-2*switchtime_MLKL),(Necr1_td-switchtime_MLKL), (Necr1_td-switchtime_MLKL/2), (Necr1_td-switchtime_MLKL/4), (Necr1_td-switchtime_MLKL/8), Necr1_td, (Necr1_td+switchtime_MLKL/8), (Necr1_td+switchtime_MLKL/4), (Necr1_td+switchtime_MLKL/2), (Necr1_td+switchtime_MLKL)], [0, 0, 0.05, 0.205, 0.340, 0.5, 0.659, 0.794, 0.95, 1], [0.025, 0.025, 0.025, 0.05, 0.075, 0.1, 0.085, 0.065, 0.05, 0.025]]).T, Necr1_obs]
    
    return ydata

"""Sloppy Cell does not have an observables component function/set. You have to create a new species and assign its in relation to other species in the model.
    Sloppy Cell ids the species as 's0', 's1'. ... The initial condiditions are 
    assigned using that notation. The observables are functions (sums) of the 
    species whose pattern contains the observable. 
    
    To get the species: Load the model, import pysb.bng. then generate_equations.
    Then, model.observables[index of the observable].species"""

#---- Initial Conditions -----
init_conditions = {}
init_conditions['Apop1'] = {'s0':600} #TNFa_0 == TNFa(brec=None) = 600 = 's0'
init_conditions['Apop2'] = {'s0':1200} #TNFa_0 == TNFa(brec=None) = 1200 = 's0'
init_conditions['Necr1'] = {'s0':1800, 's9':0, 's21':9.6e6}
#TNFa_0 == TNFa(brec=None)  = 1800 = 's0'
#FADD_0 == FADD(bDD=None, bDED1=None, bDED2=None) = 0 = 's9'
#zVAD_0 == zVad(bC8=None) = 9.6e6 = 's21'

#---- Observables ------
obs = {}
obs['Apop1'] = ('Obs_CytoC', 's789 + s792 + s793')
obs['Apop2'] = ('Obs_cPARP', 's744')
obs['Necr1'] = ('Obs_MLKL' , 's316 + s490')

#Apop1 = Experiment('Apop1')
#Apop2 = Experiment('Apop2')
Necr1 = Experiment('Necr1')

# ------Making the Data Sloppy Compatible------
expt_names = ['Necr1']
ydata = ydata_fn()
data = {}

for expt_name in expt_names:
    y_ystd = zip(ydata[expt_name][0][:,1],ydata[expt_name][0][:,2])
    t_y_std = zip(ydata[expt_name][0][:,0], y_ystd)

    data[expt_name] = {}
    data[expt_name][obs[expt_name][0]] = dict(t_y_std)

print data

Necr1.set_data(data)
