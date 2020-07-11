# pylint: disable=C0114
# pylint: disable=C0116
# pylint: disable=R0124
# pylint: disable=R0201
# pylint: disable=C0121
"Test exponential convergence"
import numpy as np
import copy
import pytest

from opensir.models import SIRX

DEFAULT_PARAMS = [
    6.2 / 8,
    1 / 8,
    0.05,
    0.05,
    5,
]  # Default parameters from Maier et al (2020)
GUANGDONG_IC = [104299972, 14, 0, 14]  # Ealing initial conditions


class TestSirx:
    """ Test class for the SIRX model """

    @pytest.fixture
    def model(self):
        return SIRX()

    def test_guangdong(self, model):
        # Check that guangdong 22 days prediction
        # is accurate with validated test case
        t = 22  # 23 days
        model.set_params(DEFAULT_PARAMS, GUANGDONG_IC)
        model.solve(t, t + 1)
        infected_end = model.fetch()[-1, 4]

        assert abs(infected_end - 2524.388983) < 0.0001

    def test_expconv(self, model):
        # Check that the numerical solution of the
        # reduction from sir-x to sir
        # converges to I ~ I_0 * exp (alpha-beta)
        # as (S,t) -> (1,0)
        # Assume small amount of infected
        n_i_0 = 1
        # Assume large population
        initial_conds = [GUANGDONG_IC[0], n_i_0, 0, n_i_0]
        # Assume very short time
        t = 1 / (24 * 3600)
        # Make kappa, kappa_0 and inf_over_test zero to reduce
        # the SIR-X to SIR
        sirx_red_params = [DEFAULT_PARAMS[0], DEFAULT_PARAMS[1], 0, 0, 0]
        # Calculate analytical and numerical solutions
        model.set_params(sirx_red_params, initial_conds)
        model.solve(t, 2)
        n_inf_num = model.fetch()[-1, 1]
        n_inf_analytical = n_inf_num * np.exp(
            (DEFAULT_PARAMS[0] - DEFAULT_PARAMS[1]) * t
        )
        perc_error = abs(n_inf_analytical - n_inf_num) / n_inf_num
        # Assert over the difference < 1e-4
        assert perc_error < 1e-3


class TestUnittestSir:
    """ Unittest class for the SIRX model """

    @pytest.fixture
    def model(self):
        return SIRX()

    def test_predict(self, model):
        t = 6
        model.set_params(DEFAULT_PARAMS, GUANGDONG_IC)
        model.solve(t, t + 1)
        m_list = model.fetch()
        # Fit the model to simulated data
        model.fit(
            m_list[:, 0],
            m_list[:, model.fit_input],
            fit_index=[False, False, True, True, True],
        )

        # Create a copy of the model
        new_model = copy.copy(model)
        # Displace the last day as the initial day in fit_attr to check
        # that predicting since t = 0 reproduces model.solve(n_days)
        new_model.fit_attr["t_obs"] = [model.fit_attr["t_obs"][0]]
        new_model.fit_attr["n_obs"] = [model.fit_attr["n_obs"][0]]
        new_list = new_model.predict(t, n_X=model.fit_attr["n_obs"][0])

        # Average absolute error in the predicted variable
        aae = abs(sum(m_list[:, model.fit_input] - new_list[:, model.fit_input])) / len(
            m_list
        )
        assert aae < 1e-9

    def test_dict_and_list_params_produce_same_results(self, model):
        t = 6
        model.set_params(DEFAULT_PARAMS, GUANGDONG_IC)
        model.solve(t, t + 1)
        sol_list = model.fetch()

        d_params = {
            "alpha": DEFAULT_PARAMS[0],
            "beta": DEFAULT_PARAMS[1],
            "kappa_0": DEFAULT_PARAMS[2],
            "kappa": DEFAULT_PARAMS[3],
            "inf_over_test": DEFAULT_PARAMS[4],
        }
        d_ic = {
            "n_S0": GUANGDONG_IC[0],
            "n_I0": GUANGDONG_IC[1],
            "n_R0": GUANGDONG_IC[2],
            "n_X0": GUANGDONG_IC[3],
        }

        model.set_params(d_params, d_ic)
        model.solve(t, t + 1)
        sol_dict = model.fetch()

        assert np.array_equal(sol_list, sol_dict)

    def test_missing_dict_param_should_raise(self, model):
        d_params = {
            "alpha": DEFAULT_PARAMS[0],
            "beta": DEFAULT_PARAMS[1],
            "kappa_0": DEFAULT_PARAMS[2],
            "kappa": DEFAULT_PARAMS[3],
            "inf_over_test": DEFAULT_PARAMS[4],
        }
        d_ic = {
            "n_S0": GUANGDONG_IC[0],
            "n_I0": GUANGDONG_IC[1],
            "n_R0": GUANGDONG_IC[2],
            "n_X0": GUANGDONG_IC[3],
        }

        bad_params = {"alpha": DEFAULT_PARAMS[0], "beta": DEFAULT_PARAMS[1]}
        # Typical error trying to initialize SIR-X with SIR parameters
        bad_ic = {
            "n_S0": GUANGDONG_IC[0],
            "n_I0": GUANGDONG_IC[1],
            "n_R0": GUANGDONG_IC[2],
        }

        model.set_params(d_params, d_ic)

        with pytest.raises(SIRX.InvalidNumberOfParametersError):
            model.set_params(d_params, bad_ic)

        with pytest.raises(SIRX.InvalidNumberOfParametersError):
            model.set_params(bad_params, d_ic)
