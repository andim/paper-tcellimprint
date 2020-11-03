import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('../custom.mplstyle')

import sys
sys.path.append('..')
from lib import *

meta = pd.read_csv('data/alpha_emerson_mincount16.csv')
meta.dropna(inplace=True)
meta['age'] = meta['Age']

def statsmodels_regression(x, y):
    x = sm.add_constant(x)
    model = sm.OLS(y,x)
    results = model.fit()
    return model, results

def plot(ax, labels, masks):
    ss_model = 0
    for label, mask in zip(labels, masks):
        x, y = meta[mask]['age'], meta[mask]['alpha']-1.0
        plot_regression(x, y, ax, logy=False, fit_slope=False, fittype='statsmodels',
                        fittransform=lambda t: (t-40)/10,
                        data_label=label, extend=10,
                        label='${1:.2f}\,(\pm {3:.2f}) {0:+.3f}\,(\pm {2:.3f}) \cdot(t-40)/10$', p_cutoff=1.0, ms=1, alpha=1.0)
        model, results = statsmodels_regression((x-40)/10.0, y)
        ss_model += np.sum((results.resid)**2)
        print(results.summary())
        
    alpha_tot = meta[mask_cmv_neg | mask_cmv_pos]['alpha']
    ss_tot = np.sum((alpha_tot - np.mean(alpha_tot))**2)
    rsq = 1-ss_model/ss_tot
    ax.text(1, 1.0, 'overall $R^2={0:.2f}$'.format(rsq),
        fontsize='x-small', transform=ax.transAxes, ha='right', va='top')


fig, axes = plt.subplots(figsize=(5.0, 2.5), ncols=2)

ax = axes[0]
mask_cmv_pos = meta['CMVpos'] == True
mask_cmv_neg = meta['CMVpos'] == False
plot(ax, ['CMV negative', 'CMV positive'], [mask_cmv_neg, mask_cmv_pos])

ax = axes[1]
mask_female = meta['Sex'] == 'Female' 
mask_male = meta['Sex'] == 'Male'
plot(ax, ['Female', 'Male'], [mask_female, mask_male])

for ax in axes:
    ax.set_ylim(0.63, 2.3)
    #lims = meta[mask]['age'].min(), meta[mask]['age'].max()
    lims = [0, 75]
    ax.set_xlim(lims)
    ax.set_xticks(lims)
    ax.set_xlabel('Age')
    ax.xaxis.set_label_coords(0.5,-0.1)
    ax.set_ylabel(r'Power-law exponent $\alpha$')
    ax.legend(loc='lower right',
          fontsize='x-small',
          bbox_to_anchor=(1.0, 1.0), bbox_transform=ax.transAxes)
fig.tight_layout()

label_axes(fig, xy=(-0.2, 1.35))
fig.savefig(figure_directory+'../figure_exponent_cmv.svg')
plt.show()
