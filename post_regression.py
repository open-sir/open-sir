"""Post regression utilities"""
import numpy as np
from sklearn.utils import resample


def ci_bootstrap(
    model, t_obs, n_I_obs, population, alpha=0.95, n_iter=1000, r0_ci=True
):
    """
    Calculates the confidence interval of the parameters
    using the random sample bootstrap method.

    Parameters
    ----------
    model : open-sir Model instance
        A Model instance, such as SIR or SIRX

    t_obs : list or np.array
        List of days where the number of infected
        where measured

    n_I_obs: list or np.array
        List with measurements of number of infected people
        in a population. Must be consistent with t_obs.

    population : integer
        Population size

    alpha : numerical scalar, optional
        Percentile of the confidence interval required

    n_iter : integer, optional
        Number of random samples that will be taken to fit the model
        and perform the bootstrapping. Use n_iter >= 1000.

    Returns
    -------
    ci : numpy.array
        Numpy array with a list of lists that contain the lower and upper
        confidence intervals of each parameter.
    p_bt : numpy.array
        list of the parameters sampled on the bootstrapping. The most
        common use of this list is to plot histograms to visualize and
        try to infer the probability density function of the parameters.

    Notes:
    -----------
    This traditional random sampling bootstrap is not a good way to bootstrap
    time-series data , baceuse the data because X(t+1) is
    correlated with X(t). In any case, it provides a reference
    case and it will can be an useful method for other types
    of models. When using this function, always compare the prediction error
    with the interval provided by the function ci_block_cv.
    """

    p0 = model.p
    w0 = model.w0

    p_bt = []
    if r0_ci:
        r0_bt = []

    # Perform bootstraping
    for i in range(0, n_iter):
        t_r, I_r = resample(t_obs, n_I_obs)
        I_r = np.array(I_r)
        # Sort and redefine arrays
        idx = np.argsort(t_r)
        t_rs = t_r[idx]
        I_rs = I_r[idx]
        # Fit the model to the sampled data
        # Update model initial conditions
        nI0_rs = I_rs[0]
        nS0_rs = population - nI0_rs
        nR0_rs = 0  # We will still assume that we don't have observations of recovered
        w0_rs = [nS0_rs, nI0_rs, nR0_rs]
        model.set_params(model.p, w0_rs)
        model.fit(t_rs, I_rs, population)
        p_bt.append(model.p)
        if r0_ci:
            r0_bt.append(model.r0)
    # Calculate percentiles value
    p_low = ((1 - alpha) / 2) * 100
    p_up = (alpha + (1 - alpha) / 2) * 100
    # From the resulting distribution, extract the
    # percentile value of the parameters

    p_bt = np.array(p_bt)

    # Construct confidence intervals
    ci = []
    for i in range(0, len(model.p)):
        ci_low = np.percentile(p_bt[:, i], p_low)
        ci_up = np.percentile(p_bt[:, i], p_up)
        ci.append([ci_low, ci_up])

    if r0_ci:
        ci.append([np.percentile(r0_bt, p_low), np.percentile(r0_bt, p_up)])

    ci = np.array(ci)
    # Reconstruct model original parameters
    model.p = p0
    model.w0 = w0

    return ci, p_bt


def ci_block_cv(model, t_obs, n_I_obs, population, lags=1, min_sample=3):
    """ Calculates the confidence interval of the model parameters
    using a block cross validation appropriate for time series
    and differential systems when the value of the states in the
    time (t+1) is not independent from the value of the states in the
    time t.

    Parameters:
    -----------

    model: a open-sir model instance
        The model needs to be initialized with the parameters and
        initial conditions

    t_obs : list or np.array
        List of days where the number of infected
        where measured

    n_I_obs: list or np.parray
        List with measurements of number of infected people
        in a population. Must be consistent with t_obs.

    population : integer
        population size

    lags : integer, optional
        Defines the number of days that will be forecasted to calculate
        the mean squared error. For example, for the prediction Xp(t) and the
        real value X(t), the mean squared error will be calculated as
        mse = 1/n_boots |Xp(t+lags)-X(t+lags)|. This provides an estimate of the
        mean deviation of the predictions after "lags" days.
        alpha: percentile of the CI required

    min_sample : integer, optional
        Number of days that will be used in the train set
        to make the first prediction.

    Returns:
    --------
    mse_avg : float
        Simple average of the mean squared error between the model
        prediction for "lags" days and the real observed value.

    mse_list : numpy array
        List of the mean squared errors using (i) points to
        predict the X(i+lags) value, with i an iterator that
        goes from n_samples+1 to the end of t_obs index.

    p_list : numpy.array
        List of the parameters sampled on the bootstrapping as a
        function of time. A common use of this list is to plot the
        mean squared error against time, to identify time periods
        where the model produces the best and worst fit to the data.
    """

    p0 = model.p
    w0 = model.w0

    # Consider at least the three first datapoints
    p_list = []
    mse_list = []  # List of mean squared errors of the prediction for the time t+1
    for i in range(min_sample - 1, len(n_I_obs) - lags):
        # Fit model to a subset of the time-series data
        model.fit(t_obs[0:i], n_I_obs[0:i], population)
        # Store the rolling parameters
        p_list.append(model.p)
        # Predict for the i+1 period
        model.solve(t_obs[i + lags], numpoints=t_obs[i + lags] + 1)
        pred = model.fetch()[:, 2]
        # Calculate mean squared errors
        mse = np.sqrt((pred[i - 1 + lags] - n_I_obs[i - 1 + lags]) ** 2)
        mse_list.append(mse)

    p_list = np.array(p_list)
    mse_list = np.array(mse_list)
    mse_avg = np.mean(mse_list)

    model.p = p0
    model.w = w0

    return mse_avg, mse_list, p_list
