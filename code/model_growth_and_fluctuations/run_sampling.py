import re
import sys
sys.path.append('..')
import numpy as np
import scipy.special
import matplotlib.pyplot as plt
import matplotlib.colors
import palettable
import pandas as pd
import glob
import os.path

from lib import *
from lib.analytical import *
from lib.fitting import *

def growthlaw(T, d, t0, gamma):
    return (-d*(T-t0) + np.log((np.exp(d * T)-1)/(np.exp(d * t0)-1))/(1 + gamma))

theta = 1e5
# C0 does not matter as it's normalized away before subsampling
C0 = 1.0
d = 0.2 
gamma = 0.1
sigma0 = 0.0
sample_size = 5e5
sigma = 0.08**.5
print('exponent', 1+(d*gamma/(1+gamma))/sigma**2)


metapath = 'data/meta.csv'
if not os.path.exists(metapath):
    Ts = np.random.uniform(0.0, 80.0, 500)
    pd.DataFrame.from_dict(dict(Age=Ts)).to_csv(metapath)
else:
    df = pd.read_csv(metapath, index_col=0)
    Ts = list(df['Age'])

for i, T in enumerate(Ts):
    outpath = 'data/ff_%g.csv.gz'%i
    if not os.path.exists(outpath):
        print(i)
        # draw initial times
        ts = np.random.uniform(low=0.0, high=T, size=int(theta*T))
        # draw initial size
        logsizes = np.log(C0)
        if sigma0 > 0:
            logsizes += np.random.normal(size=len(ts))*sigma0 - sigma0/2
        # calculate deterministic dynamics
        logsizes += growthlaw(T, d, ts, gamma)
        # draw fluctuating growth
        variance = 2*sigma**2*(T-ts)
        logsizes += np.random.normal(size=len(ts))*variance**.5 -variance/2
        sizes = np.exp(logsizes)
        # subsample
        sizes_sub = np.random.poisson(lam=sample_size * sizes/np.sum(sizes))
        mask = sizes_sub>0
        sizes_sub = sizes_sub[mask]
        ts_sub = ts[mask]
        pd.DataFrame.from_dict(dict(Age=ts_sub, counts=sizes_sub)).to_csv(outpath)

