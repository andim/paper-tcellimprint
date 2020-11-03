import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.special
import scipy.optimize
from scipy.stats import poisson, nbinom

import sys
sys.path.append('..')
from lib.fitting import *
from lib import *

nsample = 30
true_alpha = 2.15
N = 1e7
M = 5e5
cutoffs = [1, 2, 5, 10, 20, 50, 100, 200]

# negative binomial parametrization
# n = 1/alpha
# p = 1/(1+alpha*mu)
def nbinom_params(mean, alpha):
    n = 1.0/alpha
    p = 1/(1+alpha*mean)
    return n, p

if __name__ == "__main__":
    for label, subsampling_func in [('poisson', lambda c, M: poisson.rvs(c/np.sum(c)*int(M))),
                        ('nbinom alpha=0.1', lambda c, M: nbinom.rvs(*nbinom_params(c/np.sum(c)*int(M), 0.01))),
                        ('nbinom alpha=1.0', lambda c, M: nbinom.rvs(*nbinom_params(c/np.sum(c)*int(M), 1.0))),
                        ('nbinom alpha=10.0', lambda c, M: nbinom.rvs(*nbinom_params(c/np.sum(c)*int(M), 10.0))),
                        ]:
        print(label)
        alphas_mean = []
        alphas_sem = []
        for cutoff in cutoffs:
            alphas = []
            for i in range(nsample):
                c = powerlaw(N, xmin=1.0, alpha=true_alpha)
                subcounts = subsampling_func(c, M) 
                alphas.append(mle_alpha_discrete(subcounts, cmin=cutoff))
            alphas_mean.append(np.mean(alphas))
            alphas_sem.append(np.std(alphas, ddof=1)/nsample**.5)
        df = pd.DataFrame(dict(cutoff=cutoffs, mean=alphas_mean, sem=alphas_sem))
        df.to_csv('data/%s.csv'%label, index=False)
