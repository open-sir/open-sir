import model
import numpy as np
from sklearn.utils import resample

def ci_bootstrap(model, t_obs, n_I_obs, population, alpha=0.95, n_iter=1000, r0_ci = True):
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

    output: list of lists with values of predictions and
    parameters

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
        R_0_bt = []

    # Perform bootstraping
    for i in range(0,n_iter):
        t_r, I_r = resample(t_obs, n_I_obs)
        I_r = np.array(I_r)
        # Sort and redefine arrays
        idx = np.argsort(t_r)
        t_rs = t_r[idx]
        I_rs = I_r[idx]
        # Fit the model to the sampled data
        model.fit(t_rs, I_rs, population)
        p_bt.append(model.p)
        if r0_ci:
            R_0_bt.append(model.R_0)
    
    # Calculate percentiles value
    p_low = ((1-alpha)/2)*100
    p_up = (alpha + (1-alpha)/2)*100
    # From the resulting distribution, extract the
    # percentile value of the parameters

    p_bt = np.array(p_bt)

    # Construct confidence intervals
    ci = []
    for i in range(0,len(model.p)):
        ci_low = np.percentile(p_bt[:,i], p_low)
        ci_up = np.percentile(p_bt[:,i], p_up)
        ci.append([ci_low,ci_up])

    if r0_ci==True:
        ci.append([np.percentile(R_0_bt, p_low),
        np.percentile(R_0_bt, p_up)])
    
    ci = np.array(ci)
    # Reconstruct model original parameters
    model.p = p0 
    model.w0 = w0

    return ci