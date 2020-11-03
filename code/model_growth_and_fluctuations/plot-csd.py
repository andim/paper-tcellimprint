import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

import sys
sys.path.append('..')
from lib import *
from lib.fitting import *

sys.path.append('../clonesizedistributions')
from common import *

plt.style.use('../custom.mplstyle')

meta = pd.read_csv('data/meta.csv')
print(meta.shape[0])
agemin, agemax = np.floor(min(meta['Age'])), np.ceil(max(meta['Age']))
print(agemin, agemax)
agebins = np.arange(agemin+1, agemax+1, 2)
meta = meta.groupby(pd.cut(meta['Age'], bins=agebins)).head(1)
print(meta.shape[0])

norm = matplotlib.colors.Normalize(vmin=agemin, vmax=agemax-1)

fig, ax = plt.subplots(figsize=(2.25, 2.25))
meta.sort_values('Age', ascending=False, inplace=True)
for sid, row in meta[(meta['Age']>=agemin) & (meta['Age']<=agemax)].iterrows():
    age = row['Age']
    print(sid, age)
    filepath = 'data/ff_' + str(sid) +'.csv.gz'
    if os.path.exists(filepath):
        df = pd.read_csv(filepath)
        counts = df['counts']
        print(mle_alpha(counts, cmin=10))
        plot_rankfrequency(counts, color=cmap(norm(age)),
                           normalize_x=True, normalize_y=False,
                           alpha=0.5, lw=0.6,
                           scalex=1.0)

garnish(ax, agemin, agemax, cmap)

fig.tight_layout()
fig.savefig(figure_directory+'csd-aging-fullmodel.svg')

plt.show()
