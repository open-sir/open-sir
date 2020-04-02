import model
import numpy as np
from sklearn.utils import resample


def ci_bootstrap(
    model, t_obs, n_I_obs, population, alpha=0.95, n_iter=1000, r0_ci=True
):
    """
    
    Calculates the confidence interval of the parameters
    using the bootstrap method
    
    inputs:
    
    model: a open-sir model instance
    t_obs: list or np.array  of days where the number of infected
    where measured
    I_obs: vector with measurements of numbers of infected
    population: population size
    alpha: percentile of the CI required

    outputs: 

    ci: list with lower and upper confidence intervals of the parameters
    p_bt: list of the parameters sampled on the bootstrapping. The most
    common use of this list is to plot histograms to visualize and
    try to infer the probability density function of the parameters.

    disclaimer:
    This traditional bootstrap is not a good way to bootstrap 
    time-series data , baceuse the data because X(t+1) is 
    correlated with X(t). In any case, it provides a reference 
    case and it will can be an useful method for other types
    of models. 
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
        I0_rs = I_rs[0] / population
        S0_rs = (population - I0_rs) / population
        R0_rs = 0  # We will still assume that we don't have observations of recovered

        model.w0 = [S0_rs, I0_rs, R0_rs]
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

    if r0_ci == True:
        ci.append([np.percentile(r0_bt, p_low), np.percentile(r0_bt, p_up)])

    ci = np.array(ci)
    # Reconstruct model original parameters
    model.p = p0
    model.w0 = w0

    return ci, p_bt


def ci_block_cv(
    model, t_obs, n_I_obs, population, lag=1, min_sample=3, alpha=0.95, r0_ci=True
):
    """ Calculates the confidence interval of the model parameters
    using a block cross validation appropriate for time series
    and differential systems when the value of the states in the
    time (t+1) is not independent from the value of the states in the
    time t.
    
    inputs:
    
    model: a open-sir model instance
    t_obs: list or np.array  of days where the number of infected
    where measured
    I_obs: vector with measurements of numbers of infected
    population: population size
    alpha: percentile of the CI required

    outputs: #

    MSE_avg:
    MSE_list:
    p_list:

    ci: list with lower and upper confidence intervals of the parameters
    p_bt: list of the parameters sampled on the bootstrapping. The most
    common use of this list is to plot histograms to visualize and
    try to infer the probability density function of the parameters.
    """

    p0 = model.p
    w0 = model.w0

    # Consider at least the three first datapoints
    p_list = []
    MSE_list = []  # List of mean squared errors of the prediction for the time t+1
    for i in range(min_sample - 1, len(n_I_obs) - lag):
        # Fit model to a subset of the time-series data
        model.fit(t_obs[0:i], n_I_obs[0:i], population)
        # Store the rolling parameters
        p_list.append(model.p)
        # Predict for the i+1 period
        sol = model.solve(t_obs[i + lag], numpoints=t_obs[i + lag] + 1)
        # Calculate mean squared errors
        MSE = np.sqrt((population * sol[i - 1 + lag, 2] - n_I_obs[i - 1 + lag]) ** 2)
        MSE_list.append(MSE)

    p_list = np.array(p_list)
    MSE_list = np.array(MSE_list)
    MSE_avg = np.mean(MSE_list)

    model.p = p0
    model.w = w0

    return MSE_avg, MSE_list, p_list
