"""Contains class and ODE system of SIR model"""
from .model import Model


def sir(w, t, p):
    # pylint: disable=unused-argument
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
