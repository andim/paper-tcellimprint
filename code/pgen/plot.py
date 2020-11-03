import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

import sys
sys.path.append('..')
from lib import *
from lib.fitting import *

plt.style.use('../custom.mplstyle')

datadir = emerson_processeddata_directory

metadata = load_metadata_emerson(filtered=True)

bins = np.array([1, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000])
missing = []
pgens = [1e-7, 1e-8, 1e-9, 1e-10]
pgens_str = ['pgen_>1e%g'%np.log10(pgen) for pgen in pgens]
fractions = {pgen_str: dict() for pgen_str in pgens_str}
for subjectid, meta in metadata.iterrows():
    print(subjectid)
    try:
        df = pd.read_csv(emerson_processeddata_directory + subjectid + '.tsv.gz')
        df = df[(df['sequenceStatus']=='In')&(~df['zeroInsertion'])]
    except (FileNotFoundError, EOFError):
        missing.append(subjectid)
        print('missing')
        continue
    for pgen, pgen_str in zip(pgens, pgens_str):
        df[pgen_str] = df['pgen'] > pgen
    dfg = df.groupby(pd.cut(df['rank'], bins-0.1))
    for pgen_str in pgens_str:
        fractions[pgen_str][subjectid] = dfg[pgen_str].mean()
print('missing files', missing)
metadata = metadata.loc[~metadata.index.isin(missing)]

for i, rank in enumerate(bins[1:]):
    for pgen_str in fractions:
        metadata['%s_%s'%(pgen_str, rank)] = [fractions[pgen_str][subjectid].iloc[i]
                                        for subjectid, row in metadata.iterrows()]

agebinsize = 10.0
agebins = np.arange(0.0, 90.0, agebinsize)
bin_ts = agebins[:-1]+agebinsize/2
bins = np.array([1, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000])
binmids = bins[1:]
print(binmids)

df_enrichments = metadata

def aggregate(df, agebins, name):
    grouped = df.groupby(pd.cut(df['Age'], bins=agebins))
    meanfreq = grouped.agg('mean')
    meanfreq = np.array([list(meanfreq[name+'%s'%rank]) for rank in bins[1:]])
    semfreq = grouped.agg('sem')
    semfreq = np.array([list(semfreq[name+'%s'%rank]) for rank in bins[1:]])
    return meanfreq, semfreq


fig, ax = plt.subplots(figsize=(3.5, 2.8))
colors = np.asarray(palettable.matplotlib.Viridis_9.mpl_colors)
marker = ["o", "v", "^", "<", ">", "1", "2", "3", "4"]


nsizes = 8
meanfreq, semfreq = aggregate(df_enrichments, agebins, 'pgen_>1e-9_')
for i in range(nsizes):
    l, = ax.plot(bin_ts, meanfreq[i, :], '-o', c=colors[i+1], label='%g'%bins[i+1])
    ax.fill_between(bin_ts,
                meanfreq[i, :]-semfreq[i, :],
               meanfreq[i, :]+semfreq[i, :], facecolor=colors[i+1], alpha=.5, edgecolor=None)
ax.set_xlabel('Age in years (binned)')
ax.set_ylim(0.0, 0.18)
ax.set_ylabel('Fraction with $P_{gen}(\sigma)>10^{-9}$')
    
ax.legend(title='Clone size rank (binned)', ncol=2, loc='lower right')

fig.tight_layout()
plt.show()
fig.savefig(figure_directory + 'enrichment_pgen.svg')
