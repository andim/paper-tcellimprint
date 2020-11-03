import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import palettable
import pandas as pd

import sys
sys.path.append('..')
from lib import *

plt.style.use('../custom.mplstyle')

names = ['oneInsertion', 'twoInsertion']

agebinsize = 10.0
agebins = np.arange(0.0, 90.0, agebinsize)
bin_ts = agebins[:-1]+agebinsize/2
bins = np.array([1, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000])
binmids = bins[1:]
print(binmids)

df_enrichments = pd.read_csv(data_directory +'emerson-enrichments.csv', index_col=0)
for name in names:
    df_enrichments[name+'500'] = (2*df_enrichments[name+'200']+3*df_enrichments[name+'500'])/5

def aggregate(df, agebins, name):
    grouped = df.groupby(pd.cut(df['Age'], bins=agebins))
    meanfreq = grouped.agg('mean')
    meanfreq = np.array([list(meanfreq[name+'%s'%rank]) for rank in bins[1:]])
    semfreq = grouped.agg('sem')
    semfreq = np.array([list(semfreq[name+'%s'%rank]) for rank in bins[1:]])
    return meanfreq, semfreq


fig, axes = plt.subplots(figsize=(7.0, 2.8), ncols=2)
colors = np.asarray(palettable.matplotlib.Viridis_9.mpl_colors)
marker = ["o", "v", "^", "<", ">", "1", "2", "3", "4"]


nsizes = 9
for j, name in enumerate(names):
    ax = axes[j]
    meanfreq, semfreq = aggregate(df_enrichments, agebins, name)
    for i in range(1, nsizes):
        l, = ax.plot(bin_ts, meanfreq[i, :], '-o', c=colors[i], label='%g'%bins[i+1])
        ax.fill_between(bin_ts,
                    meanfreq[i, :]-semfreq[i, :],
                   meanfreq[i, :]+semfreq[i, :], facecolor=colors[i], alpha=.5, edgecolor=None)
    ax.set_xlabel('Age in years (binned)')
    ax.set_ylim(0.0, 0.068)

axes[0].set_ylabel('One insertion clones')
axes[1].set_ylabel('Two insertions clones')
axes[0].legend(title='Clone size rank (binned)', ncol=2, loc='upper right')

fig.tight_layout()
label_axes(axes)
plt.show()
fig.savefig(figure_directory+'figure_fewinsertions.svg')
