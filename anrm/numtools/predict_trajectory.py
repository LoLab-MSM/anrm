# Fits ANRM 1.0 (Irvin et. al 2013) against single-cell measurements
# of caspase reporters.

import pickle
import bayessb
import random as ra 
import numpy as np
import calibratortools as ct
import simulator_1_0 as sim
import bayes_mcmc as bmc
import matplotlib.pyplot as plt

from anrm.irvin_mod_v5_tester import model

#-----------Previously Calibrated Parameters------------
initial_position = pickle.load(open('CompII_Hypthesis_123_newtopology_2run_v40_Position.pkl'))

#----Describing Published Data-----
def ydata_fn():
    """return and array of synthesized experimental data. The array is loosely based on published experiments"""
    Apop1_td = 6.0 #six hours
    Apop2_td = 4.0 #four hours
    Necr1_td = 4.0 #four hours
    
    switchtime_CytoC = 1.0 # [hrs]
    switchtime_cPARP = 0.5 #one hour
    switchtime_MLKL = 1.0 # [hrs]
    
    Apop1_obs = ['Obs_CytoC'] #Zhang et al. Monitored CytoC (Obs_CytoC) but CytoC may not have switch behavior.
    Apop2_obs = ['Obs_cPARP']
    Necr1_obs = ['Obs_MLKL']
    
    ydata = {}
    ydata['Apop1'] = [np.array([[(Apop1_td-2*switchtime_CytoC),(Apop1_td-switchtime_CytoC), (Apop1_td-switchtime_CytoC/2), (Apop1_td-switchtime_CytoC/4), (Apop1_td-switchtime_CytoC/8), Apop1_td, (Apop1_td+switchtime_CytoC/8), (Apop1_td+switchtime_CytoC/4), (Apop1_td+switchtime_CytoC/2), (Apop1_td+switchtime_CytoC)], [0, 0, 0.05, 0.205, 0.340, 0.5, 0.659, 0.794, 0.95, 1],[0.025, 0.025, 0.025, 0.05, 0.075, 0.1, 0.085, 0.065, 0.05, 0.025]]).T, Apop1_obs]
    ydata['Apop2'] = [np.array([[(Apop2_td-2*switchtime_cPARP),(Apop2_td-switchtime_cPARP), (Apop2_td-switchtime_cPARP/2), (Apop2_td-switchtime_cPARP/4), (Apop2_td-switchtime_cPARP/8), Apop2_td, (Apop2_td+switchtime_cPARP/8), (Apop2_td+switchtime_cPARP/4), (Apop2_td+switchtime_cPARP/2), (Apop2_td+switchtime_cPARP)], [0, 0, 0.05, 0.205, 0.340, 0.5, 0.659, 0.794, 0.95, 1], [0.025, 0.025, 0.025, 0.05, 0.075, 0.1, 0.085, 0.065, 0.05, 0.025]]).T, Apop2_obs]
    ydata['Necr1'] = [np.array([[(Necr1_td-2*switchtime_MLKL),(Necr1_td-switchtime_MLKL), (Necr1_td-switchtime_MLKL/2), (Necr1_td-switchtime_MLKL/4), (Necr1_td-switchtime_MLKL/8), Necr1_td, (Necr1_td+switchtime_MLKL/8), (Necr1_td+switchtime_MLKL/4), (Necr1_td+switchtime_MLKL/2), (Necr1_td+switchtime_MLKL)], [0, 0, 0.05, 0.205, 0.340, 0.5, 0.659, 0.794, 0.95, 1], [0.025, 0.025, 0.025, 0.05, 0.075, 0.1, 0.085, 0.065, 0.05, 0.025]]).T, Necr1_obs]
    
    return ydata

#----Objective Function----------
def objective_fn(position):
    """return the value of the objective function"""
    objective = []
    for k in conditions.keys():
        ysim = solve.simulate(position, observables=True, initial_conc=conditions[k])
        PARP_MLKL_signals   = ct.extract_records(ysim, ['Obs_cPARP', 'Obs_MLKL'])
        
        if (k == 'BidKO'):
            if max(PARP_MLKL_signals[0]>0):
                td_PARP = ct.calculate_time_delay(PARP_MLKL_signals[:,0], sims.tspan)
                td_MLKL = ct.calculate_time_delay(PARP_MLKL_signals[:,1], sims.tspan)
                if td_PARP < td_MLKL:
                    objective.append(abs(td_PARP - td_MLKL))
        
        else:
            ysim_array = ct.extract_records(ysim, ydata_norm[k][1])
            ysim_norm  = ct.normalize(ysim_array, option = 1)
            ysim_tp    = ct.cubic_spline(solve.options.tspan, ysim_norm, ydata_norm[k][0][:,0]*3600)
            
            if (k == 'Necr1'):
                objective.append(np.sum((ydata_norm[k][0][:,1] - ysim_tp) ** 2 / (2 * ydata_norm[k][0][:,2])))
            
            else:
                td_PARP = ct.calculate_time_delay(PARP_MLKL_signals[:,0], sims.tspan)
                td_MLKL = ct.calculate_time_delay(PARP_MLKL_signals[:,1], sims.tspan)
                if td_MLKL < td_PARP:
                    objective.append(np.sum((ydata_norm[k][0][:,1] - ysim_tp) ** 2 / (2 * ydata_norm[k][0][:,2]))+abs(td_PARP - td_MLKL))
                else:
                    objective.append(np.sum((ydata_norm[k][0][:,1] - ysim_tp) ** 2 / (2 * ydata_norm[k][0][:,2])))
    
    return np.sum(objective)

def calculate_time_delay(signal):
    if np.isnan(np.sum(signal)):
        return None
    else:
        norm_signal = ct.normalize(signal, option = 0)
        norm_signal = norm_signal.tolist()
        idx         = norm_signal.index(min(norm_signal, key = lambda x: abs(x-0.5)))
        return ct.cubic_spline(norm_signal[idx-3:idx+3], solve.options.tspan[idx-3:idx+3], [0.5], degree = 1)

def prior(mcmc, position):
    """Distance to original parameter values"""
    
    return np.sum((position - prior_ln_mean) ** 2 / ( 2 * prior_var))


def step(mcmc):
    """Print out some statistics every 20 steps"""
    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  prior=%g  post=%g' % \
            (mcmc.iter, mcmc.sig_value, mcmc.T, float(mcmc.acceptance)/(mcmc.iter+1),
             mcmc.accept_likelihood, mcmc.accept_prior, mcmc.accept_posterior)

def prediction(mcmc, n, species_idx, pred_k, readout, ynorm, plot_samples=False):
    plt.figure()
    positions = mcmc.positions[-n:]
    accepts = mcmc.accepts[-n:]
    accept_positions = positions[accepts]
    tspan = sims.tspan
    ysamples = np.empty((len(accept_positions), len(tspan)))
    for i, pos in enumerate(accept_positions):
        ysim = solve.simulate(pos, observables = True, initial_conc=conditions[pred_k])
        ysamples[i] = ct.extract_records(ysim, readout).T
    ysamples = ysamples/ysamples.max()
    #ysamples[i] = ysim[:, species_idx] / ysim[:, species_idx].max()
    ymean = np.mean(ysamples, 0)
    ystd = np.std(ysamples, 0)
    
    if plot_samples:
        for y in ysamples:
            plt.plot(tspan, y, c='gray', alpha=.01)

    plt.plot(tspan, ymean, 'b:', linewidth=2)
    std_interval = ystd[:, None] * [+1, -1]
    plt.plot(tspan, ymean[:, None] + std_interval * 0.842, 'g-.', linewidth=2)
    plt.plot(tspan, ymean[:, None] + std_interval * 1.645, 'k-.', linewidth=2)
    plt.errorbar(ynorm[0][:,0]*3600, ynorm[0][:,1], yerr = ynorm[0][:,2], fmt = None, ecolor = 'blue')
    plt.xlim(tspan[0] - 1, tspan[-1] + 1)

#----Normalize--------------
ydata = ydata_fn()
ydata_norm = ydata.copy()
normalize = ct.normalize_array
for k in ydata_norm.keys():
    ydata_norm[k] = [normalize(ydata_norm[k][0], option = 1), ydata_norm[k][1]]

#----Initial Protein Concetrations----
conditions = {}
init_conc = {'Apop1':{'TNFa_0': 600}, 'Apop2':{'TNFa_0': 1200}, 'Necr1':{'TNFa_0':1800, 'zVad_0':9.6e6, 'FADD_0':0}, 'BidKO':{'Bid_0': 0}}
ic_params  = model.parameters_initial_conditions()
for k in init_conc.keys():
    conditions[k] = ct.initial_conditions(init_conc[k].keys(), init_conc[k].values(), ic_params)

#----Simulator Settings----
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,36000,1000) #10hrs converted to seconds (1000 timepoints)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-5
sims.atol = 1e-5

solve = sim.Solver(sims)
solve.run()

#----Bayesian and MCMC Options----
opts = bmc.MCMCOpts()
opts.nsteps = 10
opts.initial_values = np.power(10, initial_position)
#opts.initial_values = solve.initial_values
opts.likelihood_fn = objective_fn
opts.prior_fn = prior
opts.step_fn = step
opts.use_hessian = False
opts.hessian_period = opts.nsteps / 10
opts.seed = 2
opts.initial_conc = conditions
opts.estimate_params = model.parameters_rules()

# values for prior calculation
prior_mean = [p.value for p in solve.options.estimate_params]
prior_ln_mean = np.log10(prior_mean)
prior_var = 6.0

mcmc = bmc.MCMC(opts)
mcmc.run()

mixed_nsteps = opts.nsteps / 2
mixed_positions = mcmc.positions[-mixed_nsteps:]
mixed_accepts = mcmc.accepts[-mixed_nsteps:]
mixed_accept_positions = mixed_positions[mixed_accepts]
marginal_mean_pos = np.mean(mixed_accept_positions, 0)

random = np.random.RandomState(opts.seed)
sigma = 0.1;
ntimes = 20;
pred_k = 'Apop2'
yout = ['Obs_cPARP']

# show prediction for C trajectory, which was not fit to
prediction(mcmc, opts.nsteps / 2, 2, pred_k, yout, ydata_norm[pred_k], plot_samples=True)
plt.show()

