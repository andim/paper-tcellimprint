from cycler import cycler

import sys
sys.path.append('..')

from lib import *
import re
import palettable

plt.style.use('../custom.mplstyle')
datadir = 'data/'

Nsample = 1e6

numbers = re.compile(r'\d+(\.\d+)?')
def match_float(string):
    match = numbers.search(string)
    value = float(match.group(0))
    return value

fig, ax = plt.subplots(figsize=(1.85, 1.8))

colors = palettable.colorbrewer.sequential.Greens_5_r.mpl_colors
ax.set_prop_cycle(cycler('color', colors))

files = sorted(glob.glob(datadir + "earlydeath*.csv.gz"), key=match_float, reverse=True)

def load(f):
    counts = np.asarray(pd.read_csv(f, squeeze=True))
    gamma = match_float(f)
    return gamma, counts

freqs_theory = np.logspace(-6, -1, 10) 
for i, f in enumerate(files):
    gamma, counts = load(f)
    if gamma < 0.1:
        continue
    handle_data, = plot_rankfrequency(counts, normalize_y=False, ax=ax, label='%g'%gamma)
    alpha = 1.0+gamma
    fmin = 1.0/Nsample
    plaw = (freqs_theory/fmin)**-alpha
    handle_theory, = ax.plot(freqs_theory, np.sum(counts>1)*plaw, ':', color=colors[i])
    if i == 0:
        legend = plt.legend([handle_data, handle_theory], ['simulation', 'theory'], ncol=1,
                            loc='upper right', bbox_to_anchor=(1.1, 0.85))
        ax.add_artist(legend)
ax.legend(title = r'$\gamma = \theta C_0/b_0$', ncol=3,
          loc='upper right',  bbox_to_anchor=(1.1, 1.1))
ax.set_ylabel('Clone size rank')
ax.set_xlabel('Normalized clone size')
ax.set_xlim(0.5/Nsample, 0.1)
ax.set_ylim(0.9, Nsample)
ax.locator_params(numticks=10, subs=(1.0,))
fig.tight_layout()
fig.savefig(figure_directory + "model_r.svg")
plt.show()
