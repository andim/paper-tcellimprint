import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('../custom.mplstyle')

import seaborn as sns

import sys
sys.path.append('..')
from lib import *

fig, ax= plt.subplots(figsize=(2.5, 2.3))
for dataset in ['emerson', 'britanova']:
    if dataset == 'britanova':
        meta = pd.read_csv('data/alpha_britanova_mincount16.csv')
    elif dataset == 'emerson':
        meta = pd.read_csv('data/alpha_emerson_mincount16.csv')
        meta.dropna(inplace=True)
        print(meta.shape)
        meta = meta[meta['Subject ID'].isin(set(pd.read_csv(data_directory+'emerson-counts.csv')['Subject ID']))]
        print(meta.shape)
        meta['age'] = meta['Age']

    mask = meta['age'] > 0
    x, y = meta[mask]['age'], meta[mask]['alpha']-1.0
    slope, intercept, rsq = plot_regression(x, y, ax, logy=False, fit_slope=False, fittype='scipy', extend=1,
            data_label=dataset.capitalize(),
            label='${0:.4f}\pm{1:.4f}\,/$year, $r^2={2:.2f}$', p_cutoff=1.0, ms=1, alpha=1.0)
ax.set_ylim(0.63, 2.2)
lims = [0, 105 if dataset == 'britanova' else 75]
ax.set_xlim(lims)
ax.set_xticks(lims)
ax.set_xlabel('Age in years')
ax.xaxis.set_label_coords(0.5,-0.1)
ax.set_ylabel('Power-law exponent\n'+r'$\alpha$')
ax.legend(loc='upper right', fontsize='x-small', bbox_to_anchor=(0.95, 1.0), bbox_transform=fig.transFigure)
fig.tight_layout(rect=(0,0,1,0.85))
fig.savefig(figure_directory+'powerlaw_exponent.svg')
plt.show()
