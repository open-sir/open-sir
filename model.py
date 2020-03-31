import numpy as np # Numerical computing
from scipy.integrate import odeint # ODE system numerical integrator

ABSERR = 1.0e-8
RELERR = 1.0e-6
SECS=7*3600*24
NUMPOINTS=250

def call_solver(func, p, w0, t, numpoints):
    """
    Internal function to wrap the solver.
    The integrating routine *odeint* requires for parameters that were previously defined:
    * func: function to be integrated.
    * y0: vector of initial conditions of the state variables.
    * t: discrete time-steps where the solution is going to be evaluated.
    * args = (): Extra arguments to pass to function. In our case, is the vector of parameters **p**
    """

    # Call the ODE solver.
    sol = odeint(func, w0, t, args=(p,), atol=ABSERR, rtol=RELERR)
    return np.insert(sol, 0, t, axis=1)


def sir(w, t, p):
    """ SIR: Simple model of disease spread
    inputs:
    w: vector of state variables [S,I,R]
    where 
        S: Fraction of the population susceptible to the infection
        I: Fraction on the population infected
        R: Fraction of the population recovered
    t: Current time
    p: vector of parameters
    
    returns:
    f: right hand side of the system of differential equations
    """
    # Unpack state variables
    S, I, R = w
    # Unpack parameters
    alpha, beta = p
    dS_dt = -alpha*S*I
    dI_dt =  alpha*S*I - beta*I
    dR_dt = beta * I
    
    f = [dS_dt, dI_dt, dR_dt]
    return f


def sirx(w, t, p):
    """ SIR-X: Dynamic outbreaks with temporally increasing
    interventions

    inputs:
    w: vector of state variables [S,I,R,X]
    where 
        S: Fraction of the population susceptible to the infection
        I: Fraction on the population infected
        R: Fraction of the population that recovered
        X: Fraction of the population that is quarantined

    t: time
    p: vector of parameters

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


class Model:
    
    def __init__(self, p, w0):
        self.setmodel()
        self.p = p
        self.w0 = w0

    def solve(self, tf_secs, numpoints):
        tspan = np.linspace(0,tf_secs,numpoints)

        a = call_solver(self.func, self.p, self.w0,
                     tspan, numpoints)
        return a

    @property
    def R_0(self):
        """ Returns reproduction number
        R_0 = alpha/beta"""
        return self.p[0]/self.p[1]

    def fit(self, t_obs, n_I_obs, population, inplace=False):
        """ Use the Levenberg-Marquardt algorithm to fit
        the parameter alpha, as beta is assumed constant

        inputs:
        t_obs: Vector of days corresponding to the observations of number of infected people
        n_I_obs: Vector of number of infected people
        population: Size of the objective population

        Returns
        """

        secs_obs = t_obs*3600*24

        def function_handle(t, alpha, beta = self.p[1], population=population):
            p = [alpha,beta]
            I_mod = call_solver(self.func, p, self.w0,
                           t, numpoints=len(t_obs))
            n_I_mod = I_mod[:,2]*population
            return n_I_mod

        # Fit alpha
        alpha_opt, pcov = curve_fit(f = function_handle,
                        xdata = secs_obs, ydata = n_I_obs, p0=self.p[0])
        p_new = [i for i in self.p]
        p_new[0] = alpha_opt[0]
        self.p=p_new
        return
        # return p_new, pcov

    pass

class SIR(Model):
    def setmodel(self):
        self.func = sir

class SIRX(Model):
    def setmodel(self):
        self.func = sirx
