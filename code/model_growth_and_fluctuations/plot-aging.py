import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import palettable
import pandas as pd
from lib import *
import glob
plt.style.use('../custom.mplstyle')

df_enrichments = pd.read_csv('data/enrichments.csv', index_col=0)
df_enrichments.rename(columns=dict(ages='Age'), inplace=True)
fig = plotting.plot_zeroinsertion_aging(df_enrichments, 'figure_mastercurve_model', alpha=1.0, misspecification_error=1.0e-3, minrank=0)
plt.show()
