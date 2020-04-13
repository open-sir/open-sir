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
    p: vector of parameters

        alpha: transmission rate
        beta: recovery / removal rate
        kappa_0: containment rate
        kappa: quarantine rate

    returns:
    right hand side of the system of differential equation
    """
    # unpack state variable
    s, i, r, x = w  # pylint: disable=W0612
    # unpack parameters
    alpha = p[0]
    beta = p[1]
    kappa_0 = p[2]
    kappa = p[3]  # pylint: disable=W0612
    # Define ODE system
    ds_dt = -alpha * s * i - kappa_0 * s
    di_dt = alpha * s * i - beta * i - kappa_0 * i - kappa * i
    dr_dt = beta * i + kappa_0 * s
    dx_dt = (kappa_0 + kappa) * i

    return [ds_dt, di_dt, dr_dt, dx_dt]


class SIRX(Model):
    """ SIRX model definition """

    CSV_ROW = ["Days", "S", "I", "R", "X"]
    NUM_PARAMS = 5
    NUM_IC = 4
    NAME = "SIRX"

    def set_params(self, p, initial_conds):
        """ Set model parameters.

        Args:
            p (list): parameters of the model [alpha, beta, kappa_0, kappa, inf_over_test].
                 All these values should be in 1/day units.
            initial_conds (list): Initial conditions (n_S0, n_I0, n_R0, n_X0), where:

                - n_S0: Total number of susceptible to the infection
                - n_I0: Total number of infected
                - n_R0: Total number of recovered
                - n_X0: Total number of quarantined

                Note: n_S0 + n_I0 + n_R0 + n_X0 = Population

                Internally, the model initial conditions are the ratios

                - S0 = n_S0/Population
                - I0 = n_I0/Population
                - R0 = n_R0/Population
                - X0 = n_X0/Population

                which is consistent with the mathematical description
                of the SIR model.

        Returns:
            SIRX: Reference to self
        """
        self._set_params(p, initial_conds)
        self.fit_input = 4  # By default, fit against containment compartment X
        return self

    def _update_ic(self):
        """Updates i_0 = (i_0/x_0)*x_0 in the context
        of parameter fitting"""
        self.w0[1] = self.p[4] * self.w0[3]

    @property
    def _model(self):
        return sirx

    @property
    def t_inf_eff(self):
        "Returns effective infectious period"
        return 1 / sum(self.p[1:4])

    @property
    def r0_eff(self):
        "Returns effective reproduction rate"
        return self.t_inf_eff * self.p[0]

    @property
    def pcl(self):
        "Returns public containment leverage"
        return self.p[2] / (self.p[2] + self.p[3])

    @property
    def q_prob(self):
        "Returns quarantine probability"
        return (self.p[2] + self.p[3]) / sum(self.p[1:4])
