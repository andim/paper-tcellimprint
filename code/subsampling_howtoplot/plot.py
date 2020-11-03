import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

import sys
sys.path.append('..')

from lib import *
from lib.fitting import *

plt.style.use('../custom.mplstyle')

N = int(1e8)
counts = powerlaw(N, alpha=2.2)
Nc = np.sum(counts)
print('%e'%Nc)
Ms = np.logspace(5, 7, 5)
subcountss = {}
for M in Ms:
    subcounts = np.random.binomial(int(M), counts/np.sum(counts))
    subcounts = subcounts[subcounts>0]
    subcountss[M] = subcounts

fig, axes = plt.subplots(figsize=(7.0, 2.0), ncols=4)
for M in Ms:
    subcounts = subcountss[M]
    print('%e'%M)
    plot_zipf(subcounts, ax=axes[0], ls='None', marker='.', ms=1.5, label='$10^{%g}$'%round(np.log10(M), 2))
    plot_rankfrequency(subcounts, normalize_y=True, normalize_x=False, ax=axes[1])
    plot_rankfrequency(subcounts, normalize_y=True, normalize_x=True, ax=axes[2])
    plot_rankfrequency(subcounts, normalize_y=False, normalize_x=True, ax=axes[3])

ax = axes[0]
ax.set_ylabel('Density')
ax.set_xlabel('Clone size')
ax.legend(title='Sample size', ncol=2)
ax = axes[1]
ax.set_ylabel('Cumulative density')
ax.set_xlabel('Clone size')
ax = axes[2]
ax.set_ylabel('Cumulative density')
ax.set_xlabel('Normalized clone size')
ax = axes[3]
ax.set_ylabel('Rank')
ax.set_xlabel('Normalized clone size')
for ax in axes:
    ax.locator_params(numticks=9)
label_axes(fig, xy=(-0.3, 0.95))
fig.tight_layout(pad=0.2)
fig.savefig(figure_directory+'subsampled.svg')
