import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import palettable
import pandas as pd

import sys
sys.path.append('..')
from lib import *

plt.style.use('../custom.mplstyle')


agebinsize = 20.0
agebins = np.arange(0.0, 90.0, agebinsize)
bin_ts = agebins[:-1]+agebinsize/2
bins = np.array([1, 1000, 10000, 100000])
binmids = bins[1:]
print(binmids)

df_enrichments = pd.read_csv(data_directory +'emerson-enrichments.csv', index_col=0,
                             true_values=['True'], false_values=['False'])
df_enrichments.dropna(subset=['CMVpos'], inplace=True)
for nzeros in range(0,3):
    zeros = '0'*nzeros
    print(zeros)
    df_enrichments['zeroInsertion1000'+zeros] = (5*df_enrichments['zeroInsertion1000'+zeros]
                + 3*df_enrichments['zeroInsertion500'+zeros]
                + 2*df_enrichments['zeroInsertion200'+zeros])/10.0

def aggregate(df):
    grouped = df.groupby(pd.cut(df['Age'], bins=agebins))
    meanfreq = grouped.agg('mean')
    meanfreq = np.array([list(meanfreq['zeroInsertion%s'%rank]) for rank in bins[1:]])
    semfreq = grouped.agg('sem')
    semfreq = np.array([list(semfreq['zeroInsertion%s'%rank]) for rank in bins[1:]])
    return meanfreq, semfreq


fig, ax = plt.subplots(figsize=(3.5, 2.5))
colors = np.asarray(palettable.cartocolors.sequential.BluGrn_3_r.mpl_colors)

nsizes = 3
for j, df in enumerate([df_enrichments[df_enrichments['CMVpos']==True],
                        df_enrichments[df_enrichments['CMVpos']==False]]):
    meanfreq, semfreq = aggregate(df)
    for i in range(0, nsizes):
        l, = ax.plot(bin_ts, meanfreq[i, :], '-o' if j == 0 else ':o', c=colors[i], label='%g'%bins[i+1] if j==0 else None)
        if i == 0:
            if j == 0:
                lpos = l
            else:
                lneg = l
        ax.fill_between(bin_ts,
                    meanfreq[i, :]-semfreq[i, :],
                   meanfreq[i, :]+semfreq[i, :], facecolor=colors[i], alpha=.5, edgecolor=None)
ax.set_xlabel('Age in years (binned)')
ax.set_xticks(agebins[1:-1])
ax.set_ylabel('Zero insertion clones')
legend_kwargs = dict(ncol=3)
legend = plt.legend(title='Clone size rank (binned)', loc='upper right', bbox_to_anchor=(1.0, 1.05), **legend_kwargs)
ax.add_artist(legend)
ax.legend([lpos, lneg], ['positive', 'negative'], title='CMV status', loc='upper right', bbox_to_anchor=(1.0, 0.83), **legend_kwargs)
ax.set_ylim(0.0, 0.09)
ax.set_yticks(np.arange(0.0, 0.09, 0.02))
ax.set_xticks(agebins[1:-1])

fig.tight_layout()
plt.show()
fig.savefig(figure_directory+'figure_zeroinsertion_cmv.svg')
