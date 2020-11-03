from cycler import cycler

import sys
sys.path.append('..')

from lib import *
import palettable

plt.style.use('../custom.mplstyle')
datadir = 'data/'

fig, ax = plt.subplots(figsize=(3.0, 2.4))

datadir = 'data/'

def plotK(files, colors=None, **kwargs):
    lines = []
    for i, f in enumerate(files):
        print(f)
        data = np.load(f)
        theta = float(data['params'][np.argwhere(data['paramnames']=='intro')[0]])
        b0 = float(data['params'][np.argwhere(data['paramnames']=='birth')[0]])
        K = float(data['params'][np.argwhere(data['paramnames']=='K')[0]])
        r = theta/b0
        counts = data['res'][0]
        counts = counts[counts>0]
        print(max(counts))
        if colors:
            kwargs.update(color=colors[i])
        l, = plot_rankfrequency(counts, normalize_x=False, normalize_y=False, ax=ax, lw=0.8, label='%g'%K, **kwargs)
        lines.append(l)
    return tuple(lines)

files = sorted(glob.glob(datadir + "K*.npz"))
exact = plotK(files, colors=palettable.colorbrewer.sequential.Greens_5_r.mpl_colors)
ax.legend(title='K', loc='upper right', ncol=2)

ax.set_ylabel("Clone size rank")
ax.set_xlabel("Clone size")
ax.set_xscale('log')
ax.set_yscale('log')
ax.locator_params(numticks=10)

fig.tight_layout(pad=0.4, w_pad=1.0)

fig.savefig(figure_directory + "../figure_saturation.svg")
plt.show()
