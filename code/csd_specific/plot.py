import glob
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from lib import *
import palettable
plt.style.use('../custom.mplstyle')

colors = palettable.colorbrewer.sequential.Greens_7_r.mpl_colors

fig, ax = plt.subplots(figsize=(2.5, 2.5))
files = glob.glob('data/sim*.npz')
files_pb = ((np.load(f)['params'][1], f) for f in files)
for i, (pb, f) in enumerate(reversed(sorted(files_pb))):
    npz = np.load(f)
    c = npz['res']
    pb = npz['params'][1]
    print(npz['paramnames'], npz['params'])
    print(c)
    plot_rankfrequency(c, normalize_y=False, normalize_x=True, ax=ax, label=pb, color=colors[i])
ax.set_ylim(1.0, 1e4)
ax.set_xlim(2e-5, 1e-1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(loc='best', title = '$p_b$')
ax.set_ylabel('Clone size rank')
ax.set_xlabel('Normalized clone size')
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
fig.tight_layout(pad=0.2)

fig.savefig(figure_directory + '../figure_specific.svg')
plt.show()
