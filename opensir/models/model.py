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


class Model:
    """ Base model definition """

    CSV_ROW = []
    NUM_PARAMS = 4
    NUM_IC = 4

    def __init__(self):
        self.sol = None
        self.p = None
        self.pop = None
        self.w0 = None
        self.pcov = None

    @property
    def _model(self):
        raise Exception()

    def _set_params(self, p, initial_conds):
        """ Set model parameters.
        Args:
            p (list): parameters of the model. The parameters units are 1/day.
            initial_conds (list): Initial conditions, in total number of individuals.
            For instance, S0 = n_S0/population, where n_S0 is the number of subjects
            who are susceptible to the disease.

        Returns:
            Model: Reference to self
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
        """ Export the output of the model in CSV format.

        Note:
            Calling this before solve() raises an exception.

        Args:
            f: file name or descriptor
            delimiter (str): delimiter of the CSV file
        """
        if self.sol is None:
            raise Exception("Missing call to solve()")

        np.savetxt(f, self.sol, header=",".join(self.__class__.CSV_ROW), delimiter=",")

    def fetch(self):
        """ Fetch the data from the model.

        Returns:
            np.array: An array with the data. The first column is the time.
        """
        return self.sol

    def solve(self, tf_days=DAYS, numpoints=NUMPOINTS):
        """ Solve using children class model.

        Args:
            tf_days (int): number of days to simulate
            numpoints (int): number of points for the simulation.

        Returns:
            Model: Reference to self
        """
        tspan = np.linspace(0, tf_days, numpoints)
        sol = call_solver(self._model, self.p, self.w0, tspan)
        # Multiply by the population
        sol[:, 1:] *= self.pop

        self.sol = sol
        return self

    @property
    def r0(self):
        """ Returns reproduction number

        Returns:
            float: r0 (alpha/beta)
        """
        return self.p[0] / self.p[1]

    def fit(self, t_obs, n_i_obs, population, fit_index=None):
        """ Use the Levenberg-Marquardt algorithm to fit
        the parameter alpha, as beta is assumed constant

        Args:
            t_obs (np.array): Vector of days corresponding to the observations
            of number of infected people
            n_i_obs (np.array): Vector of number of infected people
            population: Size of the objective population

        Returns:
            Model: Reference to self
        """

        # if no par_index is provided, fit only the first parameter
        if fit_index is None:
            fit_index = [False for i in range(len(self.p))]
            fit_index[0] = True

        # Initial values of the parameters to be fitted
        fit_params0 = np.array(self.p)[fit_index]
        # Define fixed parameters: this set of parameters won't be fitted
        # fixed_params = self.p[fix_index]

        def function_handle(t, *par_fit, pop=population):
            params = np.array(self.p)
            params[fit_index] = par_fit
            self.p = params
            i_mod = call_solver(self._model, self.p, self.w0, t)
            return i_mod[:, 2] * pop

        # Fit parameters
        par_opt, pcov = curve_fit(
            f=function_handle, xdata=t_obs, ydata=n_i_obs, p0=fit_params0
        )
        self.p[fit_index] = par_opt
        self.pcov = pcov  # This also flags that the model was fitted
        return self
        # return p_new, pcov
