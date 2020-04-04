""" Model implementation """
import numpy as np  # Numerical computing
from scipy.integrate import odeint  # ODE system numerical integrator
from scipy.optimize import curve_fit

ABSERR = 1.0e-8
RELERR = 1.0e-6
DAYS = 7
NUMPOINTS = DAYS


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

    CSV_ROW = []
    NUM_PARAMS = 4
    NUM_IC = 4
    FUNC = None

    def __init__(self):
        self.sol = None
        self.p = None
        self.pop = None
        self.w0 = None
        self.pcov = None

    def _set_params(self, p, initial_conds):
        """ Set model parameters.
        input:
        p: parameters of the model. The parameters units are 1/day.
        initial_conds: Initial conditions, in total number of individuals.
        For instance, S0 = n_S0/population, where n_S0 is the number of subjects
        who are susceptible to the disease.
        """

        num_params = self.__class__.NUM_PARAMS
        num_ic = self.__class__.NUM_IC

        if len(p) != num_params or len(initial_conds) != num_ic:
            raise Exception(
                "Invalid number of parameters \
                             or initial conditions"
            )
        self.p = p
        self.pop = np.sum(initial_conds)
        self.w0 = initial_conds / self.pop
        return self

    def export(self, f, delimiter=","):
        """ Export the output of the model in CSV format
        Calling this before solve() raises an exception.

        input:
        f: file name or descriptor
        delimiter: delimiter of the CSV file
        """
        if self.sol is None:
            raise Exception("Missing call to solve()")

        np.savetxt(f, self.sol, header=",".join(self.__class__.CSV_ROW), delimiter=",")

    def fetch(self):
        """ Fetch the data from the model.
        The first row is the time in days
        """
        return self.sol

    def solve(self, tf_days=DAYS, numpoints=NUMPOINTS):
        """ Solve using children class model.
        input:
        tf_days: number of days to simulate
        numpoints: number of points for the simulation.

        output:
        Reference to self
        """
        tspan = np.linspace(0, tf_days, numpoints)
        sol = call_solver(self.__class__.FUNC, self.p, self.w0, tspan)
        # Multiply by the population
        sol[:, 1:] *= self.pop

        self.sol = sol
        return self

    @property
    def r0(self):
        """ Returns reproduction number
        r0 = alpha/beta"""
        return self.p[0] / self.p[1]

    def fit(self, t_obs, n_i_obs, population, fit_index=None):
        """ Use the Levenberg-Marquardt algorithm to fit
        the parameter alpha, as beta is assumed constant

        inputs:
        t_obs: Vector of days corresponding to the observations of number of infected people
        n_i_obs: Vector of number of infected people
        population: Size of the objective population

        Return
        """

        # if no par_index is provided, fit only the first parameter
        if fit_index is None:
            fit_index = [False for i in range(len(self.p))]
            fit_index[0] = True

        # Initial values of the parameters to be fitted
        fit_params0 = np.array(self.p)[fit_index]
        # Define fixed parameters: this set of parameters won't be fitted
        # fixed_params = self.p[fix_index]

        def function_handle(t, *par_fit, population=population):
            params = np.array(self.p)
            params[fit_index] = par_fit
            self.p = params
            i_mod = call_solver(self.__class__.FUNC, self.p, self.w0, t)
            return i_mod[:, 2] * population

        # Fit parameters
        par_opt, pcov = curve_fit(
            f=function_handle, xdata=t_obs, ydata=n_i_obs, p0=fit_params0
        )
        self.p[fit_index] = par_opt
        self.pcov = pcov
        return self
        # return p_new, pcov


class SIR(Model):
    """ SIR model definition """

    CSV_ROW = ["Days", "S", "I", "R"]
    NUM_PARAMS = 2
    NUM_IC = 3
    FUNC = sir

    def set_params(self, p, initial_conds):
        """ Set model parameters.
        input:
        p: parameters of the model [alpha, beta]. All these
           values should be in 1/day units.
        initial_conds: Initial conditions (n_S0, n_I0, n_R0), where:
          n_S0: Total number of susceptible to the infection
          n_I0: Toral number of infected
          n_R0: Total number of recovered
          Note n_S0 + n_I0 + n_R0 = Population

          Internally, the model initial conditions are the ratios
          S0 = n_S0/Population
          I0 = n_I0/Population
          R0 = n_R0/Population
          which is consistent with the mathematical description
          of the SIR model.

        output:
        reference to self
        """
        self._set_params(p, initial_conds)
        return self


class SIRX(Model):
    """ SIRX model definition """

    CSV_ROW = ["Days", "S", "I", "R", "X"]
    NUM_PARAMS = 4
    NUM_IC = 4
    FUNC = sirx

    def set_params(self, p, initial_conds):
        """ Set model parameters.
        input:
        p: parameters of the model [alpha, beta, kappa_0, kappa]. All these
           values should be in 1/day units.
        initial_conds: Initial conditions (S0, I0, R0, X0), where:
          n_S0: Total number of susceptible to the infection
          n_I0: Total number of infected
          n_R0: Total number of recovered
          n_X0: Total number of quarantined
          Note: n_S0 + n_I0 + n_R0 + n_X0 = Population

        Internally, the model initial conditions are the ratios
          S0 = n_S0/Population
          I0 = n_I0/Population
          R0 = n_R0/Population
          X0 = n_X0/Population
          which is consistent with the mathematical description
          of the SIR model.

        output:
        reference to self
        """
        self._set_params(p, initial_conds)
        return self
