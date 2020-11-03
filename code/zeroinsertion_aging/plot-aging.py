import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import palettable
import pandas as pd
from lib import *


plt.style.use('../custom.mplstyle')

df_enrichments = pd.read_csv(data_directory +'emerson-enrichments.csv', index_col=0)
print(df_enrichments)
fig = plotting.plot_zeroinsertion_aging(df_enrichments, 'figure_mastercurve', alpha=1.0,
        misspecification_error=1e-3, minrank=0,
        maxrank_fitting=8)
plt.show()
