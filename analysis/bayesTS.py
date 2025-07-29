from VLM.bzFRInGE import *
from VLM.bzFRInGE.FRINGEBase import ExpBase
from BZ import bbTS, bbPlot, bbGIS
# from mintpy.utils import readfile, utils as ut

from rpy2 import rinterface
from rpy2.robjects.packages import importr


log = logging.getLogger('BZ')

PATH_RES = op.join(op.expanduser('~'), 'Desktop', 'VLM')

def cgpt_regress(loc):
    import pymc3 as pm
    import numpy as np
    import matplotlib.pyplot as plt
    df = pd.read_csv(f'{PATH_RES}/{loc}_dd.csv')
    x = df['ts']
    y = df['decyr']


    # Define the changepoint model
    with pm.Model() as model:
        # Define the prior for the changepoint index
        cp_index = pm.DiscreteUniform('cp_index', lower=0, upper=len(x))

        # Define the regression coefficients and priors
        coeff_before = pm.Normal('coeff_before', mu=0, sd=1)
        coeff_after = pm.Normal('coeff_after', mu=0, sd=1)

        # Define the regression model before and after the changepoint
        regression_before = coeff_before * x[:cp_index]
        regression_after = coeff_after * x[cp_index:]

        # Define the likelihood
        sigma = pm.HalfCauchy('sigma', beta=1)
        likelihood = pm.Normal('y', mu=pm.math.concatenate([regression_before, regression_after]), sd=sigma, observed=y)

        # Perform the Bayesian inference
        trace = pm.sample(2000, tune=1000)

    # Visualize the results
    pm.plot_posterior(trace, var_names=['coeff_before', 'coeff_after', 'cp_index'], color='skyblue')




if __name__ == '__main__':
    utils = importr('utils')
    utils.install_packages('rjags')

    # Load the rjags package
    jags = importr('rjags')

