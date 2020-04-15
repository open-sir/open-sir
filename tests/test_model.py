# pylint: disable=no-self-use
"Test exponential convergence"
import pytest
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
