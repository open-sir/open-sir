import numpy as np # Numerical computing
from scipy.integrate import odeint # ODE system numerical integrator

ABSERR = 1.0e-8
RELERR = 1.0e-6
DAYS=7
NUMPOINTS=250

def sirx(w, t, p):
    """ SIR-X: Dynamic outbreaks with temporally increasing
    interventions
    
    inputs:
    w: vector of state variables [S,I,R,X]
    S: Fraction of the population susceptible to the infection
    I: Fraction on the population infected
    R: Fraction of the population that recovered
    X: Fraction of the population that is quarantined
    
    returns:
    f: right hand side of the system of differential equations
    """
    # Unpack state variables
    S, I, R, X = w
    # Unpack parameters
    alpha, beta, kappa_0, kappa = p
    dS_dt = -alpha*S*I - kappa_0*S
    dI_dt =  alpha*S*I - beta*I - kappa_0*I - kappa*I
    dR_dt = kappa_0*S + beta * I
    dX_dt = (kappa_0 + kappa) * I
    
    
    f = [dS_dt, dI_dt, dR_dt, dX_dt]
    return f

def solve_sirx(alpha, beta, kappa_0, kappa,
        S0, I0, R0, X0, days=DAYS, numpoints=NUMPOINTS):
    p = [alpha, beta, kappa_0, kappa]
    w0 = [S0, I0, R0, X0]
    secs = days * 3600 * 24
    t = [secs * float(i) / (numpoints - 1) for i in range(numpoints)]

    # Call the ODE solver.
    return odeint(sirx, w0, t, args=(p,),
              atol=ABSERR, rtol=RELERR)
