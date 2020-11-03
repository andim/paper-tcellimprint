from cycler import cycler

import sys
sys.path.append('..')

from lib import *
import re
import palettable

plt.style.use('../custom.mplstyle')
datadir = 'data/'

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

fig, ax = plt.subplots(figsize=(2.0, 1.75))

colors = palettable.colorbrewer.sequential.Blues_5_r.mpl_colors
ax.set_prop_cycle(cycler('color', colors))
files = sorted(glob.glob(datadir + "death*.npz"), key=numericalSort)

counts = np.arange(1, 1000) 
gamma = 1.0/10.0
a = np.log(1+gamma)
neutral = 1/counts * np.exp(-counts*a)
neutral /= np.sum(neutral)
neutral = 1-neutral.cumsum()
neutral, = ax.step(counts, neutral, color='k')

for f in files:
    print(f)
    data = np.load(f)
    age = float(data['params'][np.argwhere(data['paramnames']=='tlim')[0]])
    counts = np.asarray(data['res'][0])
    counts = counts[counts>0]
    plot_rankfrequency(counts, normalize_y=True, normalize_x=False, ax=ax, label='%g year'%age + ('s' if age > 1 else ''))

ax.set_xlabel("Clone size")
ax.set_ylabel("Clone size rank\n(normalized)")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1.0, 2e6)
ax.set_ylim(5e-7, 1)


ax.locator_params(numticks=10, subs=(1.0,))
legend_kwargs = dict(ncol=1)
legend = plt.legend(title='Simulation', loc='upper right', bbox_to_anchor=(1.1, 1.16), **legend_kwargs)
ax.add_artist(legend)
ax.legend([neutral], ['steady\nstate'], title='Theory', loc='upper right', bbox_to_anchor=(1.1, 0.7), **legend_kwargs)


fig.tight_layout()
fig.savefig(figure_directory + "model_steadystate.svg")
plt.show()
