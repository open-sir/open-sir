"""Contains class and ODE system of SIR model"""
import numpy as np
import matplotlib.pyplot as plt
import copy
from .model import Model, _validate_params

SIR_NUM_PARAMS = 2
SIR_NUM_IC = 3


def sir(w, t, p):
    # pylint: disable=unused-argument
    """ SIR: Simple model of disease spread
    inputs:
    w: vector of state variables [S,I,R]
    where
        S: Fraction of the population susceptible to the infection
        I: Fraction on the population infected
        R: Fraction of the population removed
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


class SIR(Model):
    """ SIR model definition """

    CSV_ROW = ["Days", "S", "I", "R"]
    PARAMS = ["alpha", "beta"]
    IC = ["n_S0", "n_I0", "n_R0"]
    NAME = "SIR"

    def predict(self, n_days=7, n_I=None, n_R=None):
        """ Predict Susceptible, Infected and Removed

        Args:
            n_days (int): number of days to predict

            n_I (int): number of infected at the last
            day of available data. If no number of
            infected is provided, the value is taken
            from the last element of the number of
            infected array on which the model was
            fitted.

            n_R (int): number of removed at the last
            day of available data. If no number of
            removed is provided, the value is set as
            the number of removed calculated by the
            SIR model as a consequence of the parameter
            fitting.

        Returns:
            np.array: Array with:
                - T: days of the predictions, where T[0] represents the last
                  day of the sample and T[1] onwards the predictions.
                - S: Predicted number of susceptible
                - I: Predicted number of infected
                - R: Predicted number of removed
        """

        # Get initial values for the predictive
        # model initial conditions
        pred_ic = self.w0 * self.pop

        if n_I is None:
            # Obtain number of infected from the last known
            # data point in the sample data
            pred_ic[1] = self.fit_attr["n_obs"][-1]
        else:
            pred_ic[1] = n_I

        if n_R is None:
            # Estimate number of recovered from the predictions
            pred_ic[2] = self.fetch()[int(self.fit_attr["t_obs"][-1]), 3]

        # Calculate new number of susceptible
        pred_ic[0] = self.pop - sum(pred_ic[1:])
        # Create a shallow copy of the model
        pred_model = copy.copy(self)
        pred_model.set_params(pred_model.p, pred_ic)

        return pred_model.solve(n_days, n_days + 1).fetch()

    def set_params(self, p, initial_conds):
        """ Set model parameters.

        Args:
            p (dict or array): parameters of the model (alpha, beta). All these
                 values should be in 1/day units. If a list is used, the order
                 of parameters is [alpha, beta].
            initial_conds (list): Initial conditions (n_S0, n_I0, n_R0), where:

                - n_S0: Total number of susceptible to the infection
                - n_I0: Toral number of infected
                - n_R0: Total number of removed

                Note n_S0 + n_I0 + n_R0 = Population

                Internally, the model initial conditions are the ratios

                - S0 = n_S0/Population
                - I0 = n_I0/Population
                - R0 = n_R0/Population

                which is consistent with the mathematical description
                of the SIR model.

                If a list is used, the order of initial conditions is
                [n_S0, n_I0, n_R0]

        Deprecated:
            This function is deprecated and will be removed soon.
            Please use :py:func:`set_parameters` and :py:func:`set_initial_conds`

        Returns:
            SIR: reference to self
        """

        super().set_params(p, initial_conds)
        self.fit_input = 2  # By default, fit against infected
        return self

    def _update_ic(self):
        """Updates initial conditions if necessary"""
        return self

    def set_parameters(self, array=None, alpha=None, beta=None):
        """ Set SIR parameters

        Args:
            array (list): list of parameters of the model ([alpha, beta])
                If set, all other arguments are ignored.
                All these values should be in 1/day units.
            alpha (float): Value of `alpha` in 1/day unit.
            beta (float): Value of `beta` in 1/day unit.

        Returns:
            SIR: Reference to self
        """
        if array:
            arr = np.array(array, dtype=float)
        else:
            arr = np.array([alpha, beta], dtype=float)

        _validate_params(arr, SIR_NUM_PARAMS)

        self.fit_input = 2  # By default, fit against infected
        self.p = arr
        return self

    def set_initial_conds(self, array=None, n_S0=None, n_I0=None, n_R0=None):
        """ Set SIR initial conditions

        Args:
            array (list): List of initial conditions [n_S0, n_I0, n_R0].
                If set, all other arguments are ignored.

                - n_S0: Total number of susceptible to the infection
                - n_I0: Total number of infected
                - n_R0: Total number of removed

                Note: n_S0 + n_I0 + n_R0 = Population

        Note:
            Internally, the model initial conditions are the ratios

            - S0 = n_S0/Population
            - I0 = n_I0/Population
            - R0 = n_R0/Population

            which is consistent with the mathematical description
            of the SIR model.

        Returns:
            SIR: Reference to self
        """
        if array:
            arr = np.array(array, dtype=float)
        else:
            arr = np.array([n_S0, n_I0, n_R0], dtype=float)

        _validate_params(arr, SIR_NUM_IC)

        self.pop = np.sum(arr)
        self.w0 = arr / self.pop
        return self

    def plot(self):
        if self.sol is None:
            raise self.InitializationError(
                "Model must be either solved or fitted before plotting"
            )

        t = self.fetch()[:, 0]
        n_i = self.fetch()[:, 2]
        n_r = self.fetch()[:, 3]
        fig, ax1 = plt.subplots(figsize=[5, 5])

        color = "tab:red"
        ax1.set_xlabel("Time / days", size=14)
        ax1.set_ylabel("Number of infected", size=14)
        ax1.plot(t, n_i, color=color, linewidth=4)
        ax1.tick_params(axis="y", labelcolor=color)
        ax1.yaxis.label.set_color(color)
        ax1.tick_params(axis="both", labelsize=13)

        ax2 = ax1.twinx()  # instantiate a second axis that shares the x coordinate

        color = "tab:blue"
        ax2.set_ylabel("Number of removed", size=14)
        ax2.plot(t, n_r, color=color, linestyle=":", linewidth=4)
        ax2.tick_params(axis="y", labelcolor=color)
        ax2.yaxis.label.set_color(color)
        ax2.tick_params(axis="both", labelsize=13)

        plt.xlim(t[0], t[-1])

        plt.title("SIR model predictions", size=16)
        plt.show()

    @property
    def _model(self):
        return sir
