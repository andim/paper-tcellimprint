import re
import sys
sys.path.append('..')
import numpy as np
import scipy.special
import matplotlib.pyplot as plt
import matplotlib.colors
import palettable
import pandas as pd
import glob
plt.style.use('../custom.mplstyle')

from lib import *
from lib.analytical import *

C = np.logspace(-3, 5, 100)
x = np.log(C)
alpha = 1.2
fig, ax = plt.subplots(figsize=(2.75, 2.0))
taus = [0.5, 1, 2, 5, 10]
norm = matplotlib.colors.Normalize(vmin=np.log(min(taus)), vmax=np.log(max(taus)))
cmap = plt.cm.get_cmap("viridis")
for i, tau in enumerate(taus):
    color = cmap(norm(np.log(tau)))
    csd = csd_alphatau(x, alpha, tau)
    ax.plot(C, csd, color=color, label=tau)
plot_referencescaling(ax, x=[10.0, 1e5], exponent=-alpha)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(min(C), max(C))
ax.set_ylim(1e-7, 2e0)
ax.set_ylabel('$P(\mathrm{Clone} > C/C_0)$')
ax.set_xlabel('Clone size $C/C_0$')
ax.locator_params(numticks=10)

ax.legend(title=r'$\tau$', loc='lower left')
fig.tight_layout(pad=0.1)
fig.savefig(figure_directory + '../figure_ffoutofsteadystate.svg')
plt.show()
