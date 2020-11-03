import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import palettable
import pandas as pd
plt.style.use('../custom.mplstyle')

import sys
sys.path.append('..')
from lib import *

agebinsize = 10.0
agebins = np.arange(0.0, 90.0, agebinsize)
ages = agebins[:-1]+agebinsize/2
bins = np.array([1, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000])

df_enrichments = pd.read_csv(data_directory+'emerson-enrichments.csv', index_col=0)
grouped = df_enrichments.groupby(pd.cut(df_enrichments['Age'], bins=agebins))
meanagg = grouped.agg('mean')
mean = np.array([list(meanagg['invdjdb%s'%rank]) for rank in bins[1:]])
semagg = grouped.agg('sem')
sem = np.array([list(semagg['invdjdb%s'%rank]) for rank in bins[1:]])
mean_nn = np.array([list(meanagg['invdjdb_lev1_%s'%rank]) for rank in bins[1:]])
sem_nn = np.array([list(semagg['invdjdb_lev1_%s'%rank]) for rank in bins[1:]])


fig, axes = plt.subplots(figsize=(5.5, 2.8), ncols=2, sharey=False)
colors = np.asarray(palettable.matplotlib.Viridis_9.mpl_colors)

nsizes = 9
for i in range(nsizes):
    l, = axes[0].plot(ages, mean[i, :], '-o', c=colors[i], label='%g'%bins[i+1])
    axes[0].fill_between(ages,
                    mean[i, :]-sem[i, :],
                    mean[i, :]+sem[i, :], facecolor=colors[i], alpha=.5, edgecolor=None)
    l_nn, = axes[1].plot(ages, mean_nn[i, :], '-o', c=colors[i], label='%g'%bins[i+1])
    axes[1].fill_between(ages,
                    mean_nn[i, :]-sem_nn[i, :],
                    mean_nn[i, :]+sem_nn[i, :], facecolor=colors[i], alpha=.5, edgecolor=None)

axes[0].set_ylabel('Fraction with exact match in VDJdb')
axes[1].set_ylabel('Fraction with close match in VDJdb')
axes[0].set_ylim(0.0, 0.05)
axes[1].set_ylim(0.0, 0.3)
for ax in axes:
    ax.set_xlabel('Age')
    ax.legend(title='Clone size rank (binned)', ncol=2)
fig.tight_layout()
label_axes(fig, xy=(-0.2, 1.0))
fig.savefig(figure_directory+'invdjdb_enrichment_ageresolved.svg')
plt.show()
