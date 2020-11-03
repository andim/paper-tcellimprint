import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('../custom.mplstyle')

import sys
sys.path.append('..')
from lib import *

agemin, agemax = 0, 110

for dataset, bins in [('britanova', np.arange(1.7, 2.7, 0.05)),
                      ('emerson', np.arange(1.7, 2.7, 0.05))
                       ]:
    if dataset == 'britanova':
        meta = pd.read_csv('data/alpha_britanova_mincount16.csv')
    else:
        meta = pd.read_csv('data/alpha_emerson_mincount16.csv')
        mask = meta['Subject ID'].isin(set(pd.read_csv(data_directory+'emerson-counts.csv')['Subject ID']))
        meta = meta[mask]
    if dataset == 'britanova':
        mask = meta['age'] > 0
        meta = meta[mask]
    y = meta['alpha']-1.0
    print(dataset, np.mean(y), np.std(y, ddof=1), str_quant(np.mean(y), np.std(y, ddof=1)/len(y)**.5),
          str_quant(np.mean(y), np.std(y, ddof=1)),
          np.median(y))
