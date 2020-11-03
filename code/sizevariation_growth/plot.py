import numpy as np
import scipy.special
import matplotlib.pyplot as plt
import palettable
plt.style.use('../custom.mplstyle')

import sys
sys.path.append('..')
from lib import *

def clonesizedistvarC0(C, C0m, gamma, sigma0):
    return (0.5 * (C/C0m)**(-2 - gamma) * np.exp(0.5 * (2 + gamma) * (1 + gamma) * sigma0**2) *
            scipy.special.erfc(((3 + 2 * gamma) * sigma0**2 - 2 * np.log(C/C0m))/(2 * np.sqrt(2) * sigma0))
            )

C = np.logspace(-3, 8, 100)
x = np.log(C)
gamma = 0.2
C0m = 1.0
fig, ax = plt.subplots(figsize=(2.75, 2.0))
sigma0s = [0.5, 1, 2]
for i, sigma0 in enumerate(sigma0s):
    pdf = clonesizedistvarC0(C, C0m=C0m, gamma=gamma, sigma0=sigma0)
    ax.plot(C, pdf, label=sigma0)
plot_referencescaling(ax, x=[1e3, 1e6], exponent=-gamma-2, factor=2e-5)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(min(C), max(C))
ax.set_ylim(1e-20, 2e0)
ax.set_ylabel('$P(C)$')
ax.set_xlabel('Clone size $C/C_0$')
ax.locator_params(numticks=10)

ax.legend(title=r'$\sigma_0$', loc='upper right')
fig.tight_layout()
fig.savefig(figure_directory+'../figure_introsizevar.svg')
plt.show()
