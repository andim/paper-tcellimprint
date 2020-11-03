import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

import sys
sys.path.append('..')
from lib.plotting import *
from lib.config import *
from lib.fitting import naive_percentage

from common import *

plt.style.use('../custom.mplstyle')

meta = pd.read_csv(data_directory + 'britanova-counts.csv')
agemin, agemax = meta[meta['age']>0]['age'].min(), meta['age'].max()
# remove cord blood samples
meta = meta[meta['age']>=agemin]
print(meta.shape)
meta.sort_values('totcount', inplace=True, ascending=False)
meta = meta.groupby('age').head(1)
print(meta.shape)


norm = matplotlib.colors.Normalize(vmin=agemin, vmax=agemax-1)

fig, axes = plt.subplots(figsize=(7.5, 2.5), ncols=3, sharex=False, sharey=False)
meta.sort_values('age', ascending=False, inplace=True)
for tup in meta[['file_name', 'age']].iterrows():
    f, age = tup[1]
    print(f, age)
    df = pd.read_csv(britanova_rawdata_directory+f+'.gz', sep='\t')
    for i, (kwargs, scalex) in enumerate([
                              (
                              dict(normalize_x=False),
                              None
                              ),
                              (
                              dict(normalize_x=True),
                              None
                              ),
                              (
                              dict(normalize_x=True),
                              lambda age: 1.0/(1.0-naive_percentage(age))
                              )
            ]):
        ax = axes[i]
        plot_rankfrequency(df['count'], ax=ax, color=cmap(norm(age)),
                           alpha=0.5, lw=0.6, 
                           scalex=scalex(age) if scalex else 1.0,
                           **kwargs
                           )

for ax in axes:
    plot_insetcolorbar(agemin, agemax, cmap, ax=ax, label='Age')
    ax.locator_params(numticks=9)
    ax.set_ylabel('Rank')

for ax in axes[1:]:
    ax.set_xlim(2e-7, 5e-1)
    ax.set_ylim(0.8, 1.5e6)

axes[0].set_xlabel('Clone size')
axes[1].set_xlabel('Normalized clone size\n(sampling depth)')
axes[2].set_xlabel('Normalized clone size\n(sampling depth & memory fraction)')

label_axes(axes, xy=(-0.2, 0.95))
fig.tight_layout()
fig.savefig(figure_directory+'csd-aging-britanova-stepbystep.svg')

plt.show()
