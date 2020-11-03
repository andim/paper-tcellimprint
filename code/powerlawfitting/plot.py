import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.special
import scipy.optimize

import sys
sys.path.append('..')
from lib.fitting import *
from lib import *
plt.style.use('../custom.mplstyle')

from run_all import true_alpha


fig, ax = plt.subplots()
for label, name in [('Poisson', 'poisson'),
                    (r'Negative Binomial, $a=1$', 'nbinom alpha=1.0'),
                    (r'Negative Binomial, $a=10$', 'nbinom alpha=10.0')]:
    df = pd.read_csv('data/%s.csv'%name)
    print(label, df)
    ax.errorbar(df['cutoff'], df['mean'], 2*np.array(df['sem']), fmt='o', label=label)
ax.axhline(true_alpha, c='k')
ax.set_xlabel('Cutoff $C_{min}$')
ax.set_ylabel(r'Power-law exponent $\alpha$')
ax.set_xscale('log')
ax.set_ylim(1.9, 3.3)
ax.legend()
fig.tight_layout()
fig.savefig(figure_directory + '../figure_powerlaw_trimming_simulations.svg')

plt.show()
