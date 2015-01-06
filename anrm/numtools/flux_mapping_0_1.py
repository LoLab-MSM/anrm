"""
    Overview
    ========
    
    Maps the chemical flux through each reaction node in a reaction network. 
    
    Reaction networks are modeled as a system of ODEs. Each ODE models the dynamics of a species. The rate change of a species is the sum of the fluxes associated with that species.  Flux (r_i) is the rate of molecule conversion by a reaction. 
    
        [Stoichiometry Matrix][Reaction_flux_vector] = d[Species_Vector]
        
        Solve for Reaction_flux_vector at each time step in the
    
    The flux_mapping program will make a graph of reaction nodes (and maybe species nodes). The edges connecting the nodes will be colored or bolded to reflect the flux associated with the nodes they connect.
    --
"""

import numpy as np
import pylab as p
import pysb.bng

from numpy.linalg import inv
from pysb.integrate import odesolve
from Enzyme_Kinetics_Example_ts import model

pysb.bng.generate_equations(model)

### Generate Stoichiometry Matrix
sp = model.species
rx = model.reactions_bidirectional

Stoichiometry_Matrix =np.zeros((len(sp), len(rx)))

for ri in range(len(rx)):
    """ Reactions with repeating reactants or product (i.e. having stoichiometry not 1) are represented by listing a species multiple times in the reactant or product tuple. These occurances add to give the correct stoichiometry."""
    for si in rx[ri]['reactants']:
        Stoichiometry_Matrix[si][ri] = Stoichiometry_Matrix[si][ri]-1
    for si in rx[ri]['products']:
        Stoichiometry_Matrix[si][ri] = Stoichiometry_Matrix[si][ri]+1

###Truncating Stoichiometry Matrix
"""There may be more species than reactions. In this case, the system is overdetermined and can be determined via Least Squares Solutions to Overdetermined Systems. 
    
    The least squares solutions to [A]x = b is found by solving ([A].T)[A]x = ([A].T)*b"""
SM_mat = np.mat(Stoichiometry_Matrix)
SM = SM_mat.T * SM_mat

###Inverse the Matrix
SMinv = inv(SM)

###First Derivative of Species vs. Time
"""FIX: We need to delete the observables from the model so that odesolve will return all species"""
t = np.linspace(0, 10, 1000)
y = odesolve(model, t)

#Foward Difference for First Derivative (2nd order accuracy)
ti = t[1]-t[0]
yi0 = np.array(list(y[0]))
yi1 = np.array(list(y[1]))
yi2 = np.array(list(y[2]))
dy = ((-3./2)*yi0 + 2*yi1 + (-1./2)*yi2)/ti

#Central Difference (2nd order accuracy)
for i in range(len(y))[1:-1]:
    yi_1 = np.array(list(y[i-1]))
    yi1 = np.array(list(y[i+1]))
    dyi = ((-1./2)*yi_1+(1./2)*yi1)/ti
    dy = np.vstack((dy, dyi))

#Reverse Difference (2nd order accuracy)
yi0 = np.array(list(y[-1]))
yi_1 = np.array(list(y[-2]))
yi_2 = np.array(list(y[-3]))
dyi = ((3./2)*yi0 - 2*yi_1 + (1./2)*yi_2)/ti

dy = np.vstack((dy, dyi))

###Find Reaction Flux....
flux = SMinv*(SM_mat.T * np.mat(dy[0]).T)

for i in range(len(dy[1:])):
    fluxi = SMinv*(SM_mat.T * np.mat(dy[i]).T)
    flux = np.hstack((flux, fluxi))

flux = np.array(flux)

###Plot
p.ion()
p.plot(t, flux[0], 'r')
#p.plot(t, flux[1])
