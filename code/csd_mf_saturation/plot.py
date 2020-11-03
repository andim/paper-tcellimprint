from cycler import cycler

import sys
sys.path.append('..')

from lib import *
import palettable

plt.style.use('../custom.mplstyle')
datadir = 'data/'

fig, ax = plt.subplots(figsize=(3.0, 2.4))


def plot(files, colors=None, **kwargs):
    lines = []
    for i, f in enumerate(files):
        print(f)
        data = np.load(f)
        theta = float(data['params'][np.argwhere(data['paramnames']=='intro')[0]])
        b0 = float(data['params'][np.argwhere(data['paramnames']=='birth')[0]])
        r = theta/b0
        counts = data['res'][0]
        counts = counts[counts>0]
        print(max(counts))
        if colors:
            kwargs.update(color=colors[i])
        l, = plot_rankfrequency(counts, normalize_x=False, normalize_y=False, ax=ax, lw=0.8, **kwargs)
        lines.append(l)
    return tuple(lines)

files = sorted(glob.glob(datadir + "exact*.npz"))
exact = plot(files, colors=palettable.colorbrewer.sequential.Greens_5_r.mpl_colors)

files = sorted(glob.glob(datadir + "mf*.npz"))
mf = plot(files, colors=palettable.colorbrewer.sequential.Blues_5_r.mpl_colors)
ax.legend([exact, mf], ['exact', 'mean field'],
           handler_map={tuple: OffsetHandlerTuple()})

ax.set_ylabel("Clone size rank")
ax.set_xlabel("Clone size")
ax.set_xscale('log')
ax.set_yscale('log')
ax.locator_params(numticks=10)

fig.tight_layout(pad=0.4, w_pad=1.0)
fig.savefig(figure_directory + "../figuremf.svg")
plt.show()
