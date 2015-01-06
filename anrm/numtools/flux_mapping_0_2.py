"""
    Overview
    ========
    
    Maps the chemical flux through each reaction node in a reaction network. 
    
    Reaction networks are modeled as a system of ODEs. Each ODE models the dynamics of a species. The rate change of a species is the sum of the fluxes associated with that species.  Flux (r_i) is the rate of molecule conversion by a reaction. 
    
        [Stoichiometry Matrix][Reaction_flux_vector] = d[Species_Vector]
        
        Solve for Reaction_flux_vector at each time step in the
    
    The flux_mapping program will make a graph of reaction nodes (and maybe species nodes). The edges connecting the nodes will be colored or bolded to reflect the flux associated with the nodes they connect.
    
    This program is our second attempt at flux mapping. The anrm contains reactions whose stoichiometries are linearly-interdependent. As a result the Stoichiometry Matrix cannot be inverted. Linear interdependence is introduced by two or more reactions who cannot be distinguished via flux balances. These reactions will be combined in this model. And the resulting flux would be the flux between species connected by a reaction or an indistinguishible set of reactions.
"""

import numpy as np
import pylab as p
import pysb.bng

from numpy.linalg import inv
from pysb.integrate import odesolve
from ANRM_ComplexI import model


def get_reaction_stoichiometry(model_reactions, reaction_index, vector_length):
    """Create a row-vector that lists the stoichiometric constants for the a particular 
        reaction. Species not involved in the reaction have stoichiometry of 0.
        
        Reactions with repeating reactants or product (i.e. having stoichiometry not 1)
        are represented by listing a species multiple times in the reactant or product
        tuple. These occurances add to give the correct stoichiometry."""
    stoich_vector = np.zeros(vector_length)
    for species_index in model_reactions[reaction_index]['reactants']:
        stoich_vector[species_index] = stoich_vector[species_index]-1
    for species_index in model_reactions[reaction_index]['products']:
        stoich_vector[species_index] = stoich_vector[species_index]+1
    return stoich_vector

pysb.bng.generate_equations(model)

### Generate Stoichiometry Matrix And matching Reaction Vector
sp = model.species
rx = model.reactions_bidirectional

t = np.linspace(0, 6000, 3000)
y = odesolve(model, t)

j=0
for i in sp:
    print j, i
    j+=1

print len(rx), " reactions"

Reaction_Vector = []
Stoichiometry_Matrix =np.array([])

#Add First Reaction Stoichiometry:
Reaction_Vector = rx[0]
Stoichiometry_Matrix = get_reaction_stoichiometry(rx, 0, len(sp))

for ri in range(len(rx))[1:]:
    #Add Next Reaction
    next_rxn = get_reaction_stoichiometry(rx, ri, len(sp))
    Stoichiometry_Matrix = np.vstack((Stoichiometry_Matrix, next_rxn))
    Reaction_Vector = np.hstack((Reaction_Vector, rx[ri]))

    # Check if the last added reaction is a stoichiometric linear combination of the other reactions. If it is, then delete it from the matrix.
    num_rows = len(Stoichiometry_Matrix)
    row_rank = np.linalg.matrix_rank(Stoichiometry_Matrix)
    if num_rows != row_rank:
        Stoichiometry_Matrix = np.delete(Stoichiometry_Matrix, np.s_[-1], 0)
        Reaction_Vector = np.delete(Reaction_Vector, np.s_[-1])
        num_rows = len(Stoichiometry_Matrix)
        row_rank = np.linalg.matrix_rank(Stoichiometry_Matrix)

print "SM.shape = ", Stoichiometry_Matrix.shape
print "Reaction_vector ", len(Reaction_Vector)

###Truncating Stoichiometry Matrix
"""There may be more species than reactions. In this case, the system is overdetermined and can be determined via Least Squares Solutions to Overdetermined Systems. 
    
    The least squares solutions to [A]x = b is found by solving ([A].T)[A]x = ([A].T)*b"""
SM_mat = np.mat(Stoichiometry_Matrix.T)
SM = SM_mat.T * SM_mat

###Inverse the Matrix
SMinv = inv(SM)

###First Derivative of Species vs. Time
"""FIX: We need to delete the observables from the model so that odesolve will return all species"""

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
for ri in range(len(Reaction_Vector)):
    p.plot(t, flux[ri], label = "r%s"%Reaction_Vector[ri]['rule'])

p.xlabel('Time [Seconds]')
p.ylabel('Flux [Molecules per second - cell]')
p.title('Reaction flux for reaction nodes in Complex I assembly')
p.legend(bbox_to_anchor = [1.3, 1.05])