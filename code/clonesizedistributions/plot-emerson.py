import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

import sys
sys.path.append('..')
from lib import *

from common import *

plt.style.use('../custom.mplstyle')

meta = pd.read_csv(data_directory + 'emerson-counts.csv')
print(meta.shape[0])
agemin, agemax = np.floor(min(meta['Age'])), np.ceil(max(meta['Age']))
print(agemin, agemax)
agebins = np.arange(agemin, agemax+1, 2)
meta.sort_values('totcount', inplace=True, ascending=False)
meta = meta.groupby(pd.cut(meta['Age'], bins=agebins)).head(1)
meta.set_index('Subject ID', inplace=True)
print(meta.shape[0])

norm = matplotlib.colors.Normalize(vmin=agemin, vmax=agemax-1)

fig, ax = plt.subplots(figsize=(2.25, 2.25))
meta.sort_values('Age', ascending=False, inplace=True)
for sid, row in meta[(meta['Age']>=agemin) & (meta['Age']<=agemax)].iterrows():
    age = row['Age']
    print(sid, age)
    filepath = emerson_processeddata_directory + sid +'.tsv.gz'
    if os.path.exists(filepath):
        df = pd.read_csv(filepath)
        mincount = df['counts'].min()
        if mincount == 1:
            plot_rankfrequency(df['counts'], color=cmap(norm(age)),
                               alpha=0.5, lw=0.6,
                               scalex=1.0/(1.0-fitting.naive_percentage(age)))

garnish(ax, agemin, agemax, cmap)

fig.tight_layout()
fig.savefig(figure_directory+'csd-aging-emerson.svg')

plt.show()
