import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

import sys
sys.path.append('..')
from lib.plotting import *
from lib.config import *

plt.style.use('../custom.mplstyle')

meta = pd.read_csv(data_directory + 'metadata-britanova.txt', sep='\t')
meta = meta[meta['age'] == 0]
print(meta.shape[0])

fig, ax = plt.subplots(figsize=(2.5, 2.5))
for tup in meta[['file_name', 'age']].iterrows():
    f, age = tup[1]
    df = pd.read_csv(britanova_rawdata_directory+f+'.gz', sep='\t')
    plot_rankfrequency(df['count'], lw=0.8)

plot_referencescaling(ax)

factor = 2.5e6
minf = 2e-7
ax.set_xlim(minf, minf*factor)
ax.set_ylim(1.0, factor)

locmaj = matplotlib.ticker.LogLocator(base=10.0, numticks=10)
ax.xaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(0.1, 1.0, 0.1))
ax.xaxis.set_minor_locator(locmin)

fig.tight_layout()
fig.savefig(figure_directory+'csd-aging-britanova-cordblood.svg')

plt.show()
