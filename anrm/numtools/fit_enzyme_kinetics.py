import pickle
import bayessb
import random as ra 
import numpy as np
import calibratortools as ct
import simulator_1_0 as sim
import bayes_mcmc as bmc
import matplotlib.pyplot as plt

from enzyme_model import model

#----Experimental Data----
"""
    ydata: dict keys = name of the experimental conditions. The items in the dict are 1x2 lists
    the first item in the list is an array of data.
    Array[:,0,:] = timepoints
    array[:,1,:] = values
    array[:,2,:] = variance
    the 3rd dimension is the observables.
    the second item in the list is the observables.
    
    init_conc: dict of experimental conditions(keys) and initial mononmer concentrations (items)
    objective_fn:
    prior:
    step:
"""

#----User Defined Functions-----
def ydata_fn():
    ydata = {}
    ydata['C'] = [np.array([[1.2e8, 2.4e8, 4.8e8, 6.0e8], [0.9e5, 1.4e5, 2.01e5, 2.0e5], [0.02e5, 0.02e5, 0.05e5, 0.2e5]])]
    return ydata

def calculate_v_0(position):
    ic_params  = model.parameters_initial_conditions()
    
    v_0 = []
    for i in ydata['C'][0][0]:
        conc_A = ct.initial_conditions(['A_0'], [i], ic_params)
        ysim = solve.simulate(position, observables = True, initial_conc=conc_A)
        Avt = ct.extract_records(ysim, ['obs_C'])
        v_0.append((1/30.0)*(-3*Avt[10]/2.0+2*Avt[11]-1*Avt[12]/2.0)[0])
     
    return v_0
                           
def objective_fn(position):
    """return the value of the objective function"""
    objective = []
    v_0 = calculate_v_0(position)
    objective.append(np.sum((ydata['C'][0][1] - v_0) ** 2 / (2 * ydata['C'][0][2])))
    return np.sum(objective)
    
def prior(mcmc, position):
    """Distance to original parameter values"""
    
    return np.sum((position - prior_ln_mean) ** 2 / ( 2 * prior_var))


def step(mcmc):
    """Print out some statistics every 20 steps"""
    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  prior=%g  post=%g' % \
            (mcmc.iter, mcmc.sig_value, mcmc.T, float(mcmc.acceptance)/(mcmc.iter+1),
             mcmc.accept_likelihood, mcmc.accept_prior, mcmc.accept_posterior)

#----Experiment Name--------
Exp_name = ('Enzyme_Kinetics_v0')

#----Data and conditions----
ydata = ydata_fn()

#----Initial Protein Concetrations----
conditions = {}
ic_params  = model.parameters_initial_conditions()

#----Simulator Settings----
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,9000,301) #seconds
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-5
sims.atol = 1e-5

solve = sim.Solver(sims)
solve.run()

#----Bayesian and MCMC Options----
opts = bmc.MCMCOpts()
opts.nsteps = 500
opts.likelihood_fn = objective_fn
opts.prior_fn = prior
opts.step_fn = step
opts.seed = ra.randint(0,1000)
#opts.initial_values = np.power(10, initial_position)
opts.initial_values = solve.initial_values
opts.initial_conc = conditions
opts.T_init = 10

# values for prior calculation
prior_mean = [p.value for p in solve.options.estimate_params]
prior_ln_mean = np.log10(prior_mean)
prior_var = 6.0

mcmc = bmc.MCMC(opts)
mcmc.run()

# print some information about the maximum-likelihood estimate parameter set
print
print '%-10s %-12s %-12s %s' % ('parameter', 'actual', 'fitted', 'log10(fit/actual)')
fitted_values = solve.cur_params(mcmc.position)[solve.estimate_idx]
for param, new_value in zip(sims.estimate_params, fitted_values):
    change = np.log10(new_value / param.value)
    values = (param.name, param.value, new_value, change)
    print '%-10s %-12.2g %-12.2g %-+6.2f' % values


# plot data
plt.ion()
initial_params = [p.value for p in sims.estimate_params]

plt.errorbar(ydata['C'][0][0], ydata['C'][0][1], ydata['C'][0][2], fmt = '%s.' % 'r', label = 'kinetics data')
yinitial = calculate_v_0(np.log10(initial_params))
plt.plot(ydata['C'][0][0], yinitial, '--b', label = 'initial')
yfinal = calculate_v_0(mcmc.position)
plt.plot(ydata['C'][0][0], yfinal, '-g', label = 'final')
plt.legend()

plt.xlabel('[A] molecules/cell')
plt.ylabel('reaction rate molecules/cell-s')
plt.title('Enzyme Kinetics')

mcmc.num_estimate = len(opts.initial_values)
hess = mcmc.calculate_hessian()

print "hessian = ", hess

plt.figure('figure 2')
from SloppyCell.ReactionNetworks import *
import scipy
#-----------Plot Eigen-spectra--------------------------
e, v = Utility.eig(hess)
e = scipy.real(e)

width = (234/len(e))**0.25 * 0.25
l = Plotting.plot_eigval_spectrum(e, offset = 0.15,
                                  widths=0.7, lw=width)




