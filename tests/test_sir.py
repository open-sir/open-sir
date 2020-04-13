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
