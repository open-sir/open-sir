# pylint: disable=C0114
# pylint: disable=C0116
# pylint: disable=R0124
# pylint: disable=R0201
# pylint: disable=C0121
"Test exponential convergence"
import numpy as np
import pytest

from opensir.models import SIR

DEFAULT_PARAMS = [0.95, 0.38]  # Default parameters from WHO
EALING_IC = [341555, 445, 0]  # Ealing initial conditions


class TestSir:
    """ Test class for the SIR model """

    @pytest.fixture
    def model(self):
        return SIR()

    def test_sixdays(self, model):
        # Check that six days prediction is accurate with validated test case
        t = 6  # 7 days
        model.set_params(DEFAULT_PARAMS, EALING_IC)
        model.solve(t, t + 1)
        removed_end = model.fetch()[-1, -1]

        assert abs(removed_end - 8360.68517) < 0.0001

    def test_expconv(self, model):
        # Check that the numerical solution converges to
        # I ~ I_0 * exp (alpha-beta) as (S,t) -> (1,0)
        # Asume small amount of infected
        n_i_0 = 1
        # Asume large population
        initial_conds = [1e8, n_i_0, 0]
        # Asume very short time
        t = 0.01
        # Calculate analytical and numerical solutions
        model.set_params(DEFAULT_PARAMS, initial_conds)
        model.solve(t, 2)
        n_inf_num = model.fetch()[-1, -1]
        n_inf_analytical = n_inf_num * np.exp(
            (DEFAULT_PARAMS[1] - DEFAULT_PARAMS[0]) * t
        )
        # Asert over the difference < 1e-4
        assert abs(n_inf_analytical - n_inf_num) < 1e-4


class TestUnittestSir:
    """ Unittest class for the SIR model """

    @pytest.fixture
    def model(self):
        return SIR()

    def test_dict_and_list_params_produce_same_results(self, model):
        t = 6
        model.set_params(DEFAULT_PARAMS, EALING_IC)
        model.solve(t, t + 1)
        m_list = model.fetch()

        d_params = {"alpha": DEFAULT_PARAMS[0], "beta": DEFAULT_PARAMS[1]}
        d_ic = {"n_S0": EALING_IC[0], "n_I0": EALING_IC[1], "n_R0": EALING_IC[2]}

        model.set_params(d_params, d_ic)
        model.solve(t, t + 1)
        m_dict = model.fetch()

        assert np.array_equal(m_list, m_dict)

    def test_missing_dict_param_should_raise(self, model):
        d_params = {"alpha": DEFAULT_PARAMS[0], "beta": DEFAULT_PARAMS[1]}
        d_ic = {"n_S0": EALING_IC[0], "n_I0": EALING_IC[1], "n_R0": EALING_IC[2]}

        bad_params = {"kappa": DEFAULT_PARAMS[0], "beta": DEFAULT_PARAMS[1]}
        bad_ic = {"n_X0": EALING_IC[0], "n_I0": EALING_IC[1], "n_R0": EALING_IC[2]}

        model.set_params(d_params, d_ic)

        with pytest.raises(SIR.InvalidNumberOfParametersError):
            model.set_params(d_params, bad_ic)

        with pytest.raises(SIR.InvalidNumberOfParametersError):
            model.set_params(bad_params, d_ic)
