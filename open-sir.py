from model import solve_sirx
import numpy as np
import matplotlib.pyplot as plt # Visualization

beta = 0.38/(3600*24) # Per day
alpha = 2.5 * beta # Reproduction number
kappa_0 = 0 # Simple SIR to begin 
kappa = 0   # Simple SIR to begin
P_Ealing = 200000 # Ealing population
I_Ealing = 200    # Infected people at 28/03/2020
S0 = (P_Ealing-I_Ealing)/P_Ealing
I0 = I_Ealing/P_Ealing
R0 = 0   # Recovered people
X0 = 0      # Quarantined people
wsol = solve_sirx([alpha, beta, kappa_0, kappa], [S0, I0, R0, X0])
# Unpack the solution
S = wsol[:,0]
I = wsol[:,1]
R = wsol[:,2]
X = wsol[:,3]
tt = np.linspace(0,7,len(S))
plt.plot(tt,S*P_Ealing)
plt.ylabel("Healthy people")
plt.xlabel("Time / d")
plt.show()
