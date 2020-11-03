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
agebins = np.arange(agemin, agemax+1, 2)
meta = meta.groupby(pd.cut(meta['age'], bins=agebins)).head(1)
print(meta.shape)


norm = matplotlib.colors.Normalize(vmin=agemin, vmax=agemax-1)

fig, ax = plt.subplots(figsize=(2.25, 2.25))
meta.sort_values('age', ascending=False, inplace=True)
for tup in meta[['file_name', 'age']].iterrows():
    f, age = tup[1]
    print(f, age)
    df = pd.read_csv(britanova_rawdata_directory+f+'.gz', sep='\t')
    plot_rankfrequency(df['count'], color=cmap(norm(age)),
                       alpha=0.5, lw=0.6,
                       scalex=1.0/(1.0-naive_percentage(age)))

garnish(ax, agemin, agemax, cmap)

fig.tight_layout()
fig.savefig(figure_directory+'csd-aging-britanova.svg')

plt.show()
