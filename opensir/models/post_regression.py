"""Post regression utilities"""
from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
from sklearn.utils import resample


@dataclass
class PredictionResults:
    """
    Class which stores model predictions,
    based on block cross-validation
    """

    pred: list  # List of out of sample predictions
    n_obs_end: float  # Last in-sample observation
    mse_avg: list  # MSE_avg after cross validation
    n_avg: list  # Sample size of each MSE_avg
    mse_seq: list  # List of lists of sequential MSE
    mse_fc: list  # List of lists of forecasting for n-days MSE
    p_cv: list  # List of parameters rolling-fitted through cross validation
    # p_bt: list  # List of parameters sampled through bootstrap

    def print_mse(self):
        """
        Prints a summary of model predictive
        performance measures for n_days forecasting
        of the fitted variable
        """

        for i, j, k in zip(self.mse_avg, self.n_avg, range(1 + self.n_avg[0])):
            print(
                "Average MSE for %.0i-day predictions = %.2f, MSE sample size = %.0i"
                % (k, i, j)
            )

    def plot_predictions(self, n_days=1):
        """
        Plot predictions and confidence intervals
        for the fitted variable
        """

        t = np.linspace(1, 1 + n_days, n_days)
        pred_low_2s = self.pred - 2 * self.mse_avg
        pred_low_s = self.pred - self.mse_avg
        pred_high_s = self.pred + self.mse_avg
        pred_high_2s = self.pred + 2 * self.mse_avg
        plt.figure(figsize=[6, 6])
        ax = plt.axes()
        ax.tick_params(axis="both", which="major", labelsize=14)

        plt.plot(t, pred_low_2s[:n_days], "b-.", linewidth=2)
        plt.plot(t, pred_low_s[:n_days], "b--", linewidth=2)
        plt.plot(t, self.pred[:n_days], "k-", linewidth=3)
        plt.plot(t, pred_high_s[:n_days], "r--", linewidth=2)
        plt.plot(t, pred_high_2s[:n_days], "r-.", linewidth=2)

        # Fancy filling
        plt.fill_between(
            t, pred_low_s[:n_days], self.pred[:n_days], alpha=0.3, color="b"
        )
        plt.fill_between(
            t, pred_low_2s[:n_days], self.pred[:n_days], alpha=0.15, color="b"
        )
        plt.fill_between(
            t, pred_high_s[:n_days], self.pred[:n_days], alpha=0.3, color="r"
        )
        plt.fill_between(
            t, pred_high_2s[:n_days], self.pred[:n_days], alpha=0.15, color="r"
        )

        plt.legend(
            [
                "-2$\sigma$ 95% IC",  # pylint: disable=W1401
                "-$\sigma$ 66% IC",  # pylint: disable=W1401
                "Prediction",
                "+$\sigma$ 66% IC",  # pylint: disable=W1401
                "+2$\sigma$ 95% IC",  # pylint: disable=W1401
            ],
            fontsize=12,
        )
        plt.plot([0, 1], [self.n_obs_end, self.pred[0]], "k-")

        plt.title("Model predictions", size=15)
        plt.xlabel("Days since the end of the sample", size=14)
        plt.ylabel("Number of infected", size=14)
        plt.show()


def _sort_resample(t_obs, n_obs):
    """
    Resamples and sort consistently an array
    of day times and cumulative number of
    infected

    Args:
        t_obs(numpy.array): time vector where the
        resampling will be performed

        n_obs(numpy.array): vector of number of
        observations, consistent with t_obs.
    """
    t_r, n_i_r = resample(t_obs, n_obs)
    n_i_r = np.array(n_i_r)
    # Sort and redefine arrays
    idx = np.argsort(t_r)
    t_rs = t_r[idx]  # Resampled sorted time
    n_i_rs = n_i_r[idx]  # Resampled sorted number of infected
    return t_rs, n_i_rs


def _percentile_to_ci(alpha, p_bt):
    """
    Input alpha and list of bootstrapped values
    for one parameter
    Output: confidence intervals of the parameters
    """

    p_low = ((1 - alpha) / 2) * 100
    p_up = (alpha + (1 - alpha) / 2) * 100
    # From the resulting distribution, extract the
    # percentile value of the parameters

    # Construct confidence intervals
    return [np.percentile(p_bt, p_low), np.percentile(p_bt, p_up)]


def rolling_avg(x_list):
    """Calculates average mean squared errors
    over block cross validation error lists"""
    err_sum = np.zeros(len(x_list))
    err_list = [[] for i in range(len(x_list))]
    n_elem = []
    for i in x_list:
        err_list.append([])
        for j, k in enumerate(i):
            err_sum[j] += k
            err_list[j].append(k)
        n_elem.append(len(i))
    avg_err = err_sum / np.linspace(len(x_list), 1, len(x_list))
    return avg_err, n_elem, err_list


class ConfidenceIntervalsMixin:
    """ Mixin with confidence interval definitions """

    def ci_bootstrap(self, alpha=0.95, n_iter=1000, r0_ci=True):
        """ Calculates the confidence interval of the parameters
        using the random sample bootstrap method.

        The model needs to be initialized and fitted prior calling
        ci_bootstrap

        Args:
            alpha (float): Percentile of the confidence interval required.

            n_iter (int): Number of random samples that will be taken to
                fit the model and perform the bootstrapping. Use n_iter >= 1000

            r0_ci (boolean): Set to True to also return the reproduction
                rate confidence interval.

        Note:
            This traditional random sampling bootstrap is not a good way to
            bootstrap time-series data , baceuse the data because X(t+1) is
            correlated with X(t). In any case, it provides a reference case and
            it will can be an useful method for other types of models. When
            using this function, always compare the prediction error with the
            interval provided by the function ci_block_cv.

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

        """

        p0 = self.p

        p_bt = []
        if r0_ci:
            r0_bt = []

        # Perform bootstraping
        for i in range(0, n_iter):  # pylint: disable=W0612
            t_rs, n_i_rs = _sort_resample(
                self.fit_attr["t_obs"], self.fit_attr["n_obs"]
            )
            w0_rs = [self.pop - n_i_rs[0], n_i_rs[0], 0]  # Still assume r0=0
            self.set_params(self.p, w0_rs)
            self.fit(t_rs, n_i_rs, self.fit_attr["fit_index"])
            p_bt.append(self.p)
            if r0_ci:
                r0_bt.append(self.r0)

        p_bt = np.array(p_bt)

        ci = []
        # Calculate and append confidence intervals for each parameters
        for i in range(len(self.p)):
            ci.append(_percentile_to_ci(alpha, p_bt[:, i]))
        # If true, calculate and append confidence interval for r0
        if r0_ci:
            ci.append(_percentile_to_ci(alpha, r0_bt))

        ci = np.array(ci)
        # Reconstruct model original parameters
        self.p = p0

        return ci, p_bt

    def block_cv(self, lags=1, min_sample=3):
        """ Calculates mean squared error of the predictions as a
        measure of model predictive performance using block
        cross validation.

        The cross-validation mean squared error can be used to
        estimate a confidence interval of model predictions.

        The model needs to be initialized and fitted
        prior calling block_cv.

        Args:
            lags (int): Defines the number of days that will be
                forecasted to calculate the mean squared error. For
                example, for the prediction Xp(t) and the real value
                X(t), the mean squared error will be calculated as mse =
                `1/n_boots |Xp(t+lags)-X(t+lags)|`. This provides an
                estimate of the mean deviation of the predictions after
                "lags" days.

            min_sample (int): Number of days that will be used in the train
                set to make the first prediction.

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
        # Consider at least the three first datapoints
        p_list = []
        mse_fc = []  # List of MSE of the prediction for the time t + i lags
        for i in range(min_sample, len(self.fit_attr["n_obs"]) + 1):
            # Fit model to a subset of the time-series data
            self.fit(
                self.fit_attr["t_obs"][0:i],
                self.fit_attr["n_obs"][0:i],
                self.fit_attr["fit_index"],
            )
            # Store the rolling parameters
            p_list.append(self.p)
            # Predict for the i + lags period
            self.solve(
                self.fit_attr["t_obs"][-1],
                numpoints=int(self.fit_attr["t_obs"][-1]) + 1,
            )
            pred = self.fetch()[:, self.fit_attr["fit_input"]]
            # Calculate mean squared errors
            if i < len(self.fit_attr["t_obs"]):
                idx = self.fit_attr["t_obs"][i - 1 + lags :].astype(int)
                mse = np.sqrt((pred[idx] - self.fit_attr["n_obs"][i - 1 + lags :]) ** 2)
                mse_fc.append(mse)

        p_list = np.array(p_list)
        mse_fc = np.array(mse_fc)
        mse_avg, n_avg, mse_seq = rolling_avg(mse_fc)

        # Generate len(t_obs) - min_sample days predictions
        self.solve(
            self.fit_attr["t_obs"][-1] + len(mse_fc),
            numpoints=int(self.fit_attr["t_obs"][-1]) + len(mse_fc) + 1,
        )

        # The first element of the prediction is out of the sample
        t_start = int(self.fit_attr["t_obs"][-1]) + 1

        # create PredictionResults dataclass
        pred_data = PredictionResults(
            pred=self.fetch()[t_start:, self.fit_attr["fit_input"]],
            n_obs_end=self.fetch()[t_start - 1, self.fit_attr["fit_input"]],
            mse_avg=mse_avg,
            n_avg=n_avg,
            mse_seq=mse_seq,
            mse_fc=mse_fc,
            p_cv=p_list,
        )

        return mse_avg, mse_seq, p_list, pred_data
