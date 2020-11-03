import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from lib import *

import matplotlib.ticker

cmap = plt.cm.get_cmap("viridis")

def garnish(ax, agemin, agemax, cmap):
    plot_insetcolorbar(agemin, agemax, cmap, ax=ax, label='Age')
    plot_referencescaling(ax)
    ax.set_xlim(2e-7, 5e-1)
    ax.set_ylim(0.8, 2.5e6)

    ax.locator_params(numticks=9)

    locmaj = matplotlib.ticker.LogLocator(base=10.0, numticks=10)
    ax.xaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(0.1, 1.0, 0.1))
    ax.xaxis.set_minor_locator(locmin)
    ax.set_xlabel('Normalized clone size')
    ax.set_ylabel('Clone size rank')
