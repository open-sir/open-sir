# pylint: disable=no-self-use
"Test exponential convergence"
import pytest
import numpy as np
from opensir.models import Model


class TestModel:
    """ Test class for the base Model """

    @pytest.fixture
    def model(self):
        """Model fixture object used by all tests"""
        return Model()

    def test_invalid_params_should_raise(self, model):
        """Tests that initializing the model with invalid params raises an
            exception.
        """
        p = [-1, 2, 3, 4]
        ic = [5, 6, 7, 8]

        with pytest.raises(Model.InvalidParameterError):
            model.set_params(p, ic)

    def test_invalid_number_params_should_raise(self, model):
        """Tests that initializing the model with invalid number of parameters
           or initial conditions raises an exception.
        """
        p = [1, 2, 3, 4]
        ic = [5, 6, 7, 8]

        with pytest.raises(Model.InvalidNumberOfParametersError):
            model.set_params([], ic)

        with pytest.raises(Model.InvalidNumberOfParametersError):
            model.set_params(p, [])

    def test_set_params_with_invalid_dict_params_should_fail(self, model):
        """Tests that calling set_params with invalid number of dict parameters
           and initial conditions raises an exception.
        """
        model.PARAMS = ["alpha", "beta"]
        model.IC = ["n_S0", "n_I0", "n_R0"]

        p = {"alpha": -1, "beta": 2}
        ic = {"n_S0": -1, "n_I0": 2, "n_R0": 3}

        with pytest.raises(Model.InvalidParameterError):
            model.set_params(p, ic)

    def test_set_params_with_invalid_number_dict_params_should_fail(self, model):
        """Tests that calling set_params with invalid number of dict parameters
           and initial conditions raises an exception.
        """
        model.PARAMS = ["alpha", "beta"]
        model.IC = ["n_S0", "n_I0", "n_R0"]

        p = {"alpha": 1, "beta": 2}
        ic = {"n_S0": 1, "n_I0": 2, "n_R0": 3}

        with pytest.raises(Model.InvalidNumberOfParametersError):
            model.set_params({}, ic)

        with pytest.raises(Model.InvalidNumberOfParametersError):
            model.set_params(p, {})

    def test_fit_parameter_consistency(self, model):
        """ Tests that running the fit function with inconsistent inputs
        between n_i_obs and t_obs, or fit_index and self.p, raise error"""

        t_obs = np.array([0, 1, 2, 3, 4])
        n_i_obs = np.array([0, 10, 20])
        # Test with a fit_index for an unrealistic
        # number of parameters
        model.p = [1, 2]
        fit_index = [True, False, False]

        # Test inconsistent dimensions of t_obs and n_i_obs
        with pytest.raises(Model.InconsistentDimensionsError):
            model.fit(t_obs, n_i_obs)
        # Testinconsistent dimnesions between fit_index and self.p
        with pytest.raises(Model.InconsistentDimensionsError):
            model.fit(t_obs[:3], n_i_obs, fit_index)

    def test_input_validity(self, model):
        """ Test that providing negative times, negative number of
        infections or non monotonic increasing times or number of
        infections outputs an error"""

        t_obs = np.array([0, 1, 2, 3])
        t_neg = np.array([-4, -3, -2, -1])
        t_dec = np.array([0, 1, 2, 1])
        t_list = [0, 1, 2, 3]

        n_obs = np.array([0, 10, 20, 35])
        n_neg = np.array([-35, 0, 10, 20])
        n_list = [0, 10, 20, 35]

        model.p = [2, 1]

        pop = 1000

        # Check for negative times
        with pytest.raises(Model.InvalidParameterError):
            model.fit(t_neg, n_obs, pop)
        # Check for non monotonically increasing time
        with pytest.raises(Model.InvalidParameterError):
            model.fit(t_dec, n_obs, pop)
        # Check for negative number of observed infections
        with pytest.raises(Model.InvalidParameterError):
            model.fit(t_obs, n_neg, pop)
        # Check that passing a list of times raises error
        with pytest.raises(Model.InvalidParameterError):
            model.fit(t_list, n_obs, pop)
        # Check that passing a list of n_obs raises an error
        with pytest.raises(Model.InvalidParameterError):
            model.fit(t_obs, n_list, pop)
