""" Model implementation """
import numpy as np  # Numerical computing
from scipy.integrate import odeint  # ODE system numerical integrator
from scipy.optimize import curve_fit

ABSERR = 1.0e-8
RELERR = 1.0e-6
DAYS = 7
NUMPOINTS = 250


def call_solver(func, p, w0, t):
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


# pylint: disable=W0613
def sir(w, t, p):
    """ SIR: Simple model of disease spread
    inputs:
    w: vector of state variables [S,I,R]
    where
        S: Fraction of the population susceptible to the infection
        I: Fraction on the population infected
        R: Fraction of the population recovered
    t: Current time
    p: vector of parameter

    returns:
    f: right hand side of the system of differential equation
    """
    # unpack state variable
    s, i, r = w  # pylint: disable=W0612
    # unpack parameter
    alpha, beta = p
    ds_dt = -alpha * s * i
    di_dt = alpha * s * i - beta * i
    dr_dt = beta * i

    return [ds_dt, di_dt, dr_dt]


def sirx(w, t, p):
    """ SIR-X: Dynamic outbreaks with temporally increasing
    intervention

    inputs:
    w: vector of state variables [S,I,R,X]
    where
        S: Fraction of the population susceptible to the infection
        I: Fraction on the population infected
        R: Fraction of the population that recovered
        X: Fraction of the population that is quarantined

    t: time
    p: vector of parameter

    returns:
    right hand side of the system of differential equation
    """
    # unpack state variable
    s, i, r, x = w  # pylint: disable=W0612
    # unpack parameter
    alpha, beta, kappa_0, kappa = p
    ds_dt = -alpha * s * i - kappa_0 * s
    di_dt = alpha * s * i - beta * i - kappa_0 * i - kappa * i
    dr_dt = kappa_0 * s + beta * i
    dx_dt = (kappa_0 + kappa) * i

    return [ds_dt, di_dt, dr_dt, dx_dt]


class Model:
    """ Base model definition """

    def __init__(self, p, w0):
        self.func = None
        self._setmodel()
        self.p = p
        self.w0 = w0

    def solve(self, tf_days=DAYS, numpoints=NUMPOINTS):
        """ Solve using children class model.
        input:
        tf_days: number of days to simulate
        numpoints: number of points for the simulation.
        """
        tspan = np.linspace(0, tf_days, numpoints)

        return call_solver(self.func, self.p, self.w0, tspan)

    # pylint: disable=R0201
    def _setmodel(self):
        raise Exception("Parent class cannot be initialized")

    @property
    def r0(self):
        """ Returns reproduction number
        r0 = alpha/beta"""
        return self.p[0] / self.p[1]

    def fit(
        self, t_obs, n_i_obs, population, fit_index=True, inplace=False,
    ):
        """ Use the Levenberg-Marquardt algorithm to fit
        the parameter alpha, as beta is assumed constant

        inputs:
        t_obs: Vector of days corresponding to the observations of number of infected people
        n_i_obs: Vector of number of infected people
        population: Size of the objective population

        Return
        """

        # if no par_index is provided, fit all model parameters
        if fit_index == True:
            fit_index = [True for i in range(len(self.p))]
        else:
            fix_index = [not i for i in fit_index]

        days_obs = t_obs

        # Initial values of the parameters to be fitted
        fit_params0 = np.array(self.p)
        # Define fixed parameters: this set of parameters won't be fitted
        # fixed_params = self.p[fix_index]

        def function_handle(t, fit_param0, fit_param1, population=population):
            params = np.array(self.p)
            fit_params = np.array([fit_param0, fit_param1])
            params[fit_index] = fit_params[fit_index]
            self.p = params
            print(fit_params)
            i_mod = call_solver(self.func, self.p, self.w0, t)
            return i_mod[:, 2] * population

        # Fit parameters
        par_opt = curve_fit(
            f=function_handle, xdata=days_obs, ydata=n_i_obs, p0=fit_params0
        )
        # p_new = np.array(self.p)
        # p_new[0] = alpha_opt[0]
        # self.p = p_new
        self.p[fit_index] = par_opt[fit_index]
        # return p_new, pcov


class SIR(Model):
    """ SIR model definition """

    def _setmodel(self):
        self.func = sir


class SIRX(Model):
    """ SIRX model definition """

    def _setmodel(self):
        self.func = sirx
