"""Post regression utilities"""
import numpy as np
from sklearn.utils import resample


def _sort_resample(time_obs, n_i_obs):
    """
    Resamples and sort consistently an array
    of day times and cumulative number of
    infected
    """
    t_r, n_i_r = resample(time_obs, n_i_obs)
    n_i_r = np.array(n_i_r)
    # Sort and redefine arrays
    idx = np.argsort(t_r)
    t_rs = t_r[idx]  # Resampled sorted time
    n_i_rs = n_i_r[idx]  # Resampled sorted number of infected
    return t_rs, n_i_rs


def _percentile_to_ci(alpha, p_bt):
    """Input alpha and list of bootstrapped values
    for one parameter
    Output: confidence intervals of the parameters """
    p_low = ((1 - alpha) / 2) * 100
    p_up = (alpha + (1 - alpha) / 2) * 100
    # From the resulting distribution, extract the
    # percentile value of the parameters

    # Construct confidence intervals
    return [np.percentile(p_bt, p_low), np.percentile(p_bt, p_up)]


class ConfidenceIntervalsMixin:
    """ Mixin with confidence interval definitions """

    def ci_bootstrap(self, t_obs, n_i_obs, population, options=None):
        """ Calculates the confidence interval of the parameters
        using the random sample bootstrap method.

        The model needs to be initialized with the parameters and
        initial conditions

        Args:
            t_obs (list or np.array): List of days where the number of infected
                where measured

            n_i_obs (list or np.array): List with measurements of number of
                infected people in a population. Must be consistent with t_obs.

            population (int): Population size

            options (dict): Random sampling bootstrapping options
                - alpha : numerical scalar
                Percentile of the confidence interval required.
                Default = 0.95
                - n_iter (int): Number of random samples that will be taken to
                fit the model and perform the bootstrapping.
                Use n_iter >= 1000. (default = 1000)
                - r0_ci (boolean): Set to True to also return the reproduction
                rate R_0 confidence interval. (default = True)

        Returns:
            tuple: tuple containing:
                - ci (numpy.array):
                    list of lists that contain the lower and upper
                    confidence intervals of each parameter.
                - p_bt (numpy.array):
                    list of the parameters sampled on the bootstrapping.
                    The most common use of this list is to plot histograms to
                    visualize and try to infer the probability density function of
                    the parameters.

        Note:
            This traditional random sampling bootstrap is not a good way to
            bootstrap time-series data , baceuse the data because X(t+1) is
            correlated with X(t). In any case, it provides a reference case and
            it will can be an useful method for other types of models. When
            using this function, always compare the prediction error with the
            interval provided by the function ci_block_cv.
        """

        # If no options provided, use default confidence interval of 95%
        if options is None:
            options = {"alpha": 0.95, "n_iter": 1000, "r0_ci": True}

        p0 = self.p

        p_bt = []
        if options["r0_ci"]:
            r0_bt = []

        # Perform bootstraping
        for i in range(0, options["n_iter"]):  # pylint: disable=W0612
            t_rs, n_i_rs = _sort_resample(t_obs, n_i_obs)
            w0_rs = [population - n_i_rs[0], n_i_rs[0], 0]  # Still assume r0=0
            self.set_params(self.p, w0_rs)
            self.fit(t_rs, n_i_rs, population)
            p_bt.append(self.p)
            if options["r0_ci"]:
                r0_bt.append(self.r0)

        p_bt = np.array(p_bt)

        ci = []
        # Calculate and append confidence intervals for each parameters
        for i in range(len(self.p)):
            ci.append(_percentile_to_ci(options["alpha"], p_bt[:, i]))
        # If true, calculate and append confidence interval for r0
        if options["r0_ci"]:
            ci.append(_percentile_to_ci(options["alpha"], r0_bt))

        ci = np.array(ci)
        # Reconstruct model original parameters
        self.p = p0

        return ci, p_bt

    def ci_block_cv(self, t_obs, n_i_obs, population, options=None):
        """ Calculates the confidence interval of the model parameters
        using a block cross validation appropriate for time series
        and differential systems when the value of the states in the
        time (t+1) is not independent from the value of the states in the
        time t.

        The model needs to be initialized with the parameters and
        initial conditions

        Args:
            t_obs (list or np.array): List of days where the number of infected
                where measured

            n_i_obs (list or np.parray): List with measurements of number of
                infected people in a population. Must be consistent with t_obs.

            population (int): population size

            options (dict): Time bootstrapping options
                - lags (int): Defines the number of days that will be
                forecasted to calculate the mean squared error. For
                example, for the prediction Xp(t) and the real value
                X(t), the mean squared error will be calculated as mse =
                `1/n_boots |Xp(t+lags)-X(t+lags)|`. This provides an
                estimate of the mean deviation of the predictions after
                "lags" days.  (default = 1)
                - min_sample (int): Number of days that will be used in the train
                set to make the first prediction. (default = 3)

        Returns:
            tuple: tuple containing:
                - mse_avg (float):
                    Simple average of the mean squared error between the model
                    prediction for "lags" days and the real observed value.
                - mse_list (numpy.array):
                    List of the mean squared errors using (i) points to
                    predict the X(i+lags) value, with i an iterator that goes from
                    n_samples+1 to the end of t_obs index.
                - p_list (numpy.array):
                    List of the parameters sampled on the bootstrapping as
                    a function of time. A common use of this list is to plot the
                    mean squared error against time, to identify time periods where
                    the model produces the best and worst fit to the data.

        """

        p0 = self.p
        w0 = self.w0

        if options is None:
            options = {"lags": 1, "min_sample": 3}

        lags = options["lags"]

        # Consider at least the three first datapoints
        p_list = []
        mse_list = []  # List of mean squared errors of the prediction for the time t+1
        for i in range(options["min_sample"] - 1, len(n_i_obs) - lags):
            # Fit model to a subset of the time-series data
            self.fit(t_obs[0:i], n_i_obs[0:i], population)
            # Store the rolling parameters
            p_list.append(self.p)
            # Predict for the i+1 period
            self.solve(t_obs[i + lags], numpoints=int(t_obs[i + lags]) + 1)
            pred = self.fetch()[:, 2]
            # Calculate mean squared errors
            mse = np.sqrt((pred[i - 1 + lags] - n_i_obs[i - 1 + lags]) ** 2)
            mse_list.append(mse)

        p_list = np.array(p_list)
        mse_list = np.array(mse_list)
        mse_avg = np.mean(mse_list)

        self.p = p0
        self.w = w0

        return mse_avg, mse_list, p_list
