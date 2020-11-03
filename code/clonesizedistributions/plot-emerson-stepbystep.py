import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

import sys
sys.path.append('..')
from lib import *
from lib.fitting import *

from common import *

plt.style.use('../custom.mplstyle')

meta = pd.read_csv(data_directory + 'emerson-counts.csv')
print(meta.shape[0])
agemin, agemax = np.floor(min(meta['Age'])), np.ceil(max(meta['Age']))
print(agemin, agemax)
agebins = np.arange(agemin, agemax+1, 1)
meta.sort_values('totcount', inplace=True, ascending=False)
meta = meta.groupby(pd.cut(meta['Age'], bins=agebins)).head(1)
meta.set_index('Subject ID', inplace=True)
print(meta.shape[0])

norm = matplotlib.colors.Normalize(vmin=agemin, vmax=agemax-1)

fig, axes = plt.subplots(figsize=(7.5, 2.5), ncols=3, sharex=False, sharey=False)
meta.sort_values('Age', ascending=False, inplace=True)
for sid, row in meta[(meta['Age']>=agemin) & (meta['Age']<=agemax)].iterrows():
    age = row['Age']
    print(sid, age)
    filepath = emerson_processeddata_directory + sid +'.tsv.gz'
    if os.path.exists(filepath):
        df = pd.read_csv(filepath)
        mincount = df['counts'].min()
        if mincount == 1:
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
                plot_rankfrequency(df['counts'], ax=ax, color=cmap(norm(age)),
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


label_axes(axes, labels=('D', 'E', 'F'), xy=(-0.2, 0.95))
fig.tight_layout()
fig.savefig(figure_directory+'csd-aging-emerson-stepbystep.svg')

plt.show()
