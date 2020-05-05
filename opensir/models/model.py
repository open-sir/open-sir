""" Model implementation """
import numpy as np  # Numerical computing
from scipy.integrate import odeint  # ODE system numerical integrator
from scipy.optimize import curve_fit

from opensir.models.post_regression import ConfidenceIntervalsMixin

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


class Model(ConfidenceIntervalsMixin):
    """ Base model definition """

    CSV_ROW = []
    PARAMS = []
    IC = []
    NAME = None

    def __init__(self):
        self.sol = None
        self.p = None
        self.pop = None
        self.w0 = None
        self.fit_input = None
        self.fit_attr = None  # Dict set after the first call of fit

    class InvalidParameterError(Exception):
        """Raised when an initial parameter of a value is not correct"""

    class InvalidNumberOfParametersError(Exception):
        """Raised when the number of initial parameters is not correct"""

    class InconsistentDimensionsError(Exception):
        """Raised when the length of the days array is
        not equal to the dimension of the observed cases, or if the length
        of fit_index has a length different than the length
        of the parameter array self.p"""

    class InitializationError(Exception):
        """Raised when a function executed violating the
        logical sequence of the Open-SIR pipeline"""

    @property
    def _model(self):
        raise Exception()

    def set_params(self, p, initial_conds):
        """ Set model parameters.
        Args:
            p (dict or list): parameters of the model. The parameters units are 1/day,
                and should be >= 0.
            initial_conds (dict or list): Initial conditions, in total number of individuals.
                For instance, S0 = n_S0/population, where n_S0 is the number of subjects
                who are susceptible to the disease.

        Note:
            Support for list of params and initial conditions will be deprecated
            in future releases.

        Returns:
            Model: Reference to self
        """

        params = self.PARAMS
        ic = self.IC

        if isinstance(p, dict):
            p = [p[d] for d in params if d in p.keys()]

        if isinstance(initial_conds, dict):
            initial_conds = [initial_conds[d] for d in ic if d in initial_conds.keys()]

        num_params = len(params)
        num_ic = len(ic)

        try:
            for param in p:
                assert param >= 0
        except:
            raise self.InvalidParameterError()

        if len(p) != num_params or len(initial_conds) != num_ic:
            raise self.InvalidNumberOfParametersError()

        self.p = p
        self.pop = np.sum(initial_conds)
        self.w0 = initial_conds / self.pop
        return self

    def export(self, f, suppress_header=False, delimiter=","):
        """ Export the output of the model in CSV format.

        Note:
            Calling this before solve() raises an exception.

        Args:
            f: file name or descriptor
            suppress_header (boolean): Set to true to suppress the CSV header
            delimiter (str): delimiter of the CSV file
        """
        if self.sol is None:
            raise Exception("Missing call to solve()")

        kwargs = {"delimiter": delimiter}

        if not suppress_header:
            kwargs["comments"] = ""
            kwargs["header"] = ",".join(self.CSV_ROW)

        np.savetxt(f, self.sol, **kwargs)

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
            float:
            .. math::
                R_0 = \\alpha/\\beta
        """
        return self.p[0] / self.p[1]

    def fit(self, t_obs, n_obs, fit_index=None):
        """ Use the Levenberg-Marquardt algorithm to fit
        model parameters consistent with True entries in the fit_index list.

        Args:
            t_obs (numpy.ndarray): Vector of days corresponding to
                the observations of number of infected people.
                Must be a non-decreasing array.
            n_obs (numpy.nparray): Vector which contains the observed
                epidemiological variable to fit the model against.
                It must be consistent with t_obs and with the initial
                conditions defined when building the model and using the
                set_parameters and set_initial_conds function.
                The model fit_input attribute defines against which
                epidemiological variable the fitting will be performed.
            fit_index (list of booleans , optional): this list must be
                of the same size of the number of  parameters of the model.
                The parameter p[i] will be fitted if fit_index[i] = True. Otherwise,
                the parameter will be fixed. By default, fit will only fit the first
                parameter of p, p[0].

        Returns:
            Model: Reference to self
        """

        # If t_obs or n_obs are not numpy arrays, raise error
        if (not isinstance(t_obs, np.ndarray)) | (not isinstance(n_obs, np.ndarray)):
            raise self.InvalidParameterError("t_obs and n_obs must be a numpy arrays")

        # Raise error if any time or number of infected is negative
        if any(t_obs < 0) | any(n_obs < 0):
            raise self.InvalidParameterError(
                "Time and number of infected must be non-negative"
            )

        # Raise error if any time or number of infected lists
        # is non monotonically increasing
        if any(np.diff(t_obs) < 0):
            raise self.InvalidParameterError(
                "The array of times must be non decreasing"
            )

        # Check consistent dimensions
        if len(t_obs) != len(n_obs):
            raise self.InconsistentDimensionsError(Exception)
        # if no par_index is provided, fit only the first parameter
        if fit_index is None:
            fit_index = [False for i in range(len(self.p))]
            fit_index[0] = True
        elif len(fit_index) != len(self.p):
            raise self.InconsistentDimensionsError(Exception)

        # Initialize self.fit_attr in the first call of fit
        if self.fit_attr is None:
            self.fit_attr = {
                "fit_input": self.fit_input,  # Is not None because it depends on the model
                "fit_index": None,
                "t_obs": None,
                "n_obs": None,
                "pcov": None,
            }

        # Check consistent inputs for fitting
        # for i in
        # Set post-fit model attributes
        self.fit_attr["fit_index"] = fit_index
        if self.fit_attr["t_obs"] is None:
            self.fit_attr["t_obs"] = t_obs
        if self.fit_attr["n_obs"] is None:
            self.fit_attr["n_obs"] = n_obs

        # Initial values of the parameters to be fitted
        fit_params0 = np.array(self.p)[self.fit_attr["fit_index"]]
        # Define fixed parameters: this set of parameters won't be fitted
        # fixed_params = self.p[fix_index]

        def function_handle(t, *par_fit, pop=self.pop):
            params = np.array(self.p)
            params[self.fit_attr["fit_index"]] = par_fit
            self.p = params
            self._update_ic()  # Updates IC if necessary. For example, i_o/x_0 for SIR-X
            sol_mod = call_solver(self._model, self.p, self.w0, t)
            return sol_mod[:, self.fit_attr["fit_input"]] * pop

        # Fit parameters
        # Ensure non-negativity and a loose upper bound
        bounds = (np.zeros(len(fit_params0)), np.ones(len(fit_params0)) * 100)
        # return p_new, pcov
        par_opt, pcov = curve_fit(
            f=function_handle, xdata=t_obs, ydata=n_obs, p0=fit_params0, bounds=bounds
        )
        self.p[self.fit_attr["fit_index"]] = par_opt
        self.fit_attr["pcov"] = pcov  # This also flags that the model was fitted
        return self

    def _update_ic(self):
        """ updates initial conditions if necessary """
        return self


def _validate_params(arr, num_params):
    if np.isnan(arr).any() or len(arr) != num_params:
        raise Model.InvalidNumberOfParametersError

    if sum(arr < 0) > 0:
        raise Model.InvalidParameterError
