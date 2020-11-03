import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('../custom.mplstyle')

import sys
sys.path.append('..')
from lib import *

agemin, agemax = 6, 74

for dataset, bins in [('britanova', np.arange(1.7, 2.7, 0.05)),
                      ('britanova_cordblood', np.arange(1.7, 2.7, 0.1)),
                      ('emerson', np.arange(1.7, 2.7, 0.05))
                       ]:
    if dataset in ['britanova_cordblood', 'britanova']:
        meta = pd.read_csv('data/alpha_britanova_mincount16.csv')
    else:
        meta = pd.read_csv('data/alpha_emerson_mincount16.csv')
        #mask = meta['Subject ID'].isin(set(pd.read_csv(data_directory+'emerson-counts.csv')['Subject ID']))
        mask = (meta['Age'] > agemin) & (meta['Age'] < agemax)
        meta = meta[mask]
    if dataset == 'britanova':
        mask = meta['age'] > 0
        mask &= (meta['age'] > agemin) & (meta['age'] < agemax)
        meta = meta[mask]
    elif dataset == 'britanova_cordblood':
        mask = meta['age'] == 0
        meta = meta[mask]
    y = meta['alpha']-1.0
    bins -= 1.0
    fig, ax = plt.subplots(figsize=(1, 0.5))
    ax.hist(y, bins=bins)
    print(dataset, np.mean(y), np.std(y, ddof=1), str_quant(np.mean(y), np.std(y, ddof=1)/len(y)**.5),
          str_quant(np.mean(y), np.std(y, ddof=1)),
          np.median(y), scipy.stats.iqr(y))
    ax.text(1.0, 0.95, r'$\alpha$', va='top', ha='right', transform=ax.transAxes)
    ax.axvline(np.mean(y), c='k')
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    fig.savefig(figure_directory + 'powerlaw_exponent_%s_hist.svg'%dataset)
