from pysb import *
from pysb.integrate import odesolve
from pylab import linspace, plot, xlabel, ylabel, show

# A simple model with a reversible binding rule

Model()

# Declare the monomers
Monomer('A'  ,  ['bf'])
Monomer('B' ,   ['bf'])
Monomer('C' ,   ['bf'])


Parameter('A_0' ,  60000) # 6000 corresponds to 100ng/ml TNFa
Parameter('B_0' ,  8000) # 4800 receptors per cell
Parameter('k1',1e-6)
Parameter('k2',1e-3)
Parameter('k3',0.1)

Initial(A(bf=None), A_0)
Initial(B(bf=None), B_0)


Rule('rxn1', A(bf=None) + B(bf=None) <> A(bf = 1)%B(bf = 1), k1, k2)
Rule('rxn2', A(bf=1)%B(bf=1) >> B(bf=None) + C(bf=None), k3)

"""
# Simulate the model through 40 seconds
time = linspace(0, 40, 100)
print "Simulating..."
x = odesolve(model, time)
# Plot the trajectory of LR
plot(time, x['LR'])
xlabel('Time (seconds)')
ylabel('Amount of LR')
show()
"""