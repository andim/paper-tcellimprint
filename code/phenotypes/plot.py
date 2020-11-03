import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import matplotlib
plt.style.use('../custom.mplstyle')

import sys
sys.path.append('..')
from lib import *
from lib.fitting import *


regression_kwargs = dict(logy=True, ms=2, extend=1, fittype='scipy',
        label='fit: ${0:.0f} \pm {1:.0f}$ years')

fig, axes = plt.subplots(figsize=(6.5, 6.0), ncols=2, nrows=2, sharex=False, sharey=False)

## Chu et al.
meta = pd.read_csv(data_directory+'metadata-chu.csv')

dfs_by_phenotype = {}
for ind, row in meta[(meta['subject'] == 1) & ((meta['date'] == '2011-09-16') |
                                               (meta['date'] == '2011-09-15'))].iterrows():
    df = pd.read_csv(chu_rawdata_directory+row['path'], sep='\t', header=None, usecols=(0, 1, 2, 27, 30),
                     skiprows=1, names=['nucleotide', 'aminoAcid', 'counts', 'n1Insertion', 'n2Insertion'])

    df['freqs'] = df['counts']/df['counts'].sum()
    dfs_by_phenotype[row['phenotype']] = df
dfm = multimerge(dfs_by_phenotype.values(), 'nucleotide', suffixes=dfs_by_phenotype.keys())

ax = axes[0, 1]
for ind in dfs_by_phenotype.keys():
    plot_rankfrequency(dfm['freqs_'+ind], normalize_y=False, label=('Unsorted' if ind == 'PBMC' else ind), normalize_x=False, ax=ax)
#plot_rankfrequency(dfm['freqs_PBMC'], normalize_y=False,
#                                    scalex=2.6, label='Unsorted shifted', ax=ax)
mask = dfm['freqs_Memory'].isna()
plot_rankfrequency(dfm[mask]['freqs_Naive'], normalize_y=False, normalize_x=False,
                        label='Naive filtered', ax=ax)

plot_referencescaling(ax=ax, factor=0.5)
ax.legend(loc='upper right')
ax.set_xlim(2e-6)
ax.set_ylim(0.9, 5e5)
ax.set_xlabel('Normalized clone size')
ax.set_ylabel('Clone size rank')

ax = axes[1,0]
fmin = 9e-7 
mask = (dfm['freqs_PBMC']>fmin) & (dfm['freqs_Memory']>fmin)
dfmm = dfm[mask]
density_scatter(dfmm['freqs_PBMC'], dfmm['freqs_Memory'],
                ax=ax, s=1, bins=20,
                trans=lambda x: np.log(x),
                rasterized=True,
                cmap='viridis')
ax.set_xscale('log')
ax.set_yscale('log')
x = np.array([fmin, 1e-2])
ax.plot(x, x, 'k-', label='y = x', lw=2.5)
ax.plot(x, x*2.6, 'g-', label='y = 2.6 x', lw=2.5)
ax.set_xlim(fmin*1.1, 5e-3)
ax.set_ylim(fmin*1.1, 5e-3)
#ax.legend(loc='upper left')
ax.set_xlabel('Fraction Unsorted')
ax.set_ylabel('Fraction Memory')

## Britanova

memoryearly = pd.read_csv(data_directory + 'memoryearly.csv', sep=',', skiprows=1)
naive = pd.read_csv(data_directory + 'naive-percentage-britanova.csv', index_col=0)
naive.dropna(inplace=True)
ax = axes[1, 1]
x = naive['Age']
y = naive['Naïve,%']/100
slope, intercept, rsq = plot_regression(x, y, ax,
                                        plot_ci=True, data_label='Britanova',
                                        **regression_kwargs)

x = memoryearly['Age group']/12.
y = 1.0-memoryearly['median']/100
mask = x<6
x = np.array(x[mask])
y = np.array(y[mask])
slope_early, intercept_early, rsq = plot_regression(x, y, ax,
            plot_ci=False, data_label='Shearer (median)', **regression_kwargs)

series = pd.Series(data=[slope, np.exp(intercept),
                         slope_early, np.exp(intercept_early)],
                   index=['slope', 'intercept',
                         'slope_early', 'intercept_early'])
print(series)
series.to_csv(data_directory+'naive_fit.csv', header=False)

ax.legend(loc='lower left')
ax.set_xlabel('Age in years')
ax.set_ylabel('Fraction naive cells')

dataset, Nsample = ('britanova', 5e5)
ax = axes[0,0]
df_counts = pd.read_csv(data_directory + dataset + '-counts.csv', index_col=1)
df_counts.rename(str.lower, axis='columns', inplace=True)
mask = (df_counts['age'] > 0) & (~df_counts['sub_singletons'].isna())

dfm = pd.merge(df_counts, naive, right_index=True, left_index=True)
dfm = dfm[dfm['age']>0]
dfm.dropna(inplace=True)

x, y = dfm['Naïve,%']/100, (dfm['sub_singletons']/Nsample)
slope, intercept, rsq = plot_regression(x, y, ax=ax, plot_ci=True)
ax.text(0.1, 0.9, '$r^2 = {0:.2f}$'.format(rsq), fontsize='x-small')
ax.set_xlabel('Fraction naive\n(flow cytometry)')
ax.set_ylabel('Fraction singletons\n(TCR sequencing)')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

##

label_axes(fig, xy=(-0.22, 0.95))

fig.tight_layout()
fig.savefig(figure_directory + '../figure_naive.svg')
plt.show()
