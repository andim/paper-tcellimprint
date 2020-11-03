import pandas as pd
import numpy as np
from functools import reduce
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from lib import *
plt.style.use('../custom.mplstyle')

for britanova in [False, True]:
    mincounts = range(1, 51)
    if britanova:
        dfs = [pd.read_csv('data/alpha_britanova_mincount%i.csv'%i, index_col=0, usecols=(2,7))
               for i in mincounts]
    else:
        dfs = [pd.read_csv('data/alpha_emerson_mincount%i.csv'%i, index_col=0, usecols=(4,5))
           for i in mincounts]
    alphas = pd.concat(dfs, axis=1)
    alphas.columns = mincounts
    print(alphas)

    if britanova:
        meta = pd.read_csv('data/alpha_britanova_mincount1.csv', index_col=0)
        # filter cord blood samples
        alphas = alphas[~alphas.index.isin(set(meta[meta['age']==0]['sample_id']))]

    prng = np.random.RandomState(1234)
    fig, ax = plt.subplots()
    for name in prng.choice(alphas.index, 50):
        ax.plot(alphas.loc[name]-1, c='k', lw=.5)
    ax.plot(alphas.mean()-1, lw=3)
    inset = ax.inset_axes((0.6, 0.6, 0.38, 0.38))
    print(alphas.columns.max())
    inset.plot(alphas.var(), label='Variance')
    inset.set_yscale('log')
    inset.set_xlabel('Cutoff $C_{min}$')
    inset.set_ylabel('Variance')
    inset.set_xlim(1.0, alphas.columns.max()-1)
    inset.set_ylim(9e-3, 4e-1)
    ax.set_xlim(1.0)
    ax.set_ylim(0.5, 3.5)
    ax.set_xlabel('Cutoff $C_{min}$')
    ax.set_ylabel(r'Power-law exponent $\alpha$')
    fig.tight_layout()
    fig.savefig(figure_directory + 'powerlaw_trimming_' + ('britanova' if britanova else 'emerson') + '.svg')
    plt.show()
