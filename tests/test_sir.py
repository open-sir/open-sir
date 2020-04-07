# pylint: disable=C0114
# pylint: disable=C0116
# pylint: disable=R0124
# pylint: disable=C0121
"Test exponential convergence"
# from models.sir import SIR
import numpy as np
from models.sir import SIR


def test_sixdays():
    # Check that six days prediction is accurate with validated test case
    p = [0.95, 0.38]  # Default parameters from WHO
    nw0 = [341555, 445, 0]  # Ealing initial conditions
    t = 6  # 7 days
    my_sir = SIR()
    my_sir.set_params(p, nw0)
    my_sir.solve(t, t + 1)
    removed_end = my_sir.fetch()[-1, -1]

    assert abs(removed_end - 8360.68517) < 0.0001


def test_expconv():
    # Check that the numerical solution converges to
    # I ~ I_0 * exp (alpha-beta) as (S,t) -> (1,0)
    # Use default parameters
    p = [0.95, 0.38]
    # Asume small amount of infected
    n_i_0 = 1
    # Asume large population
    initial_conds = [1e8, n_i_0, 0]
    # Asume very short time
    t = 0.01
    # Calculate analytical and numerical solutions
    my_sir = SIR()
    my_sir.set_params(p, initial_conds)
    my_sir.solve(t, 2)
    n_inf_num = my_sir.fetch()[-1, -1]
    n_inf_analytical = n_inf_num * np.exp((p[1] - p[0]) * t)
    # Asert over the difference < 1e-4
    assert abs(n_inf_analytical - n_inf_num) < 1e-4
