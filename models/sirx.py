"""Contains class and ODE system of SIR-X model"""
from .model import Model


def sirx(w, t, p):
    # pylint: disable=unused-argument
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


class SIRX(Model):
    """ SIRX model definition """

    CSV_ROW = ["Days", "S", "I", "R", "X"]
    NUM_PARAMS = 4
    NUM_IC = 4

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

    @property
    def _model(self):
        return sirx
