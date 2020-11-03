import numpy as np
import pandas as pd
import sys
sys.path.append('..')
from lib.fitting import *
from lib.config import *
from lib import *

mincount = int(sys.argv[1])

def alpha(path):
    print(path)
    try:
        counts = pd.read_csv(path)['counts']
        alpha = mle_alpha_discrete(counts, cmin=mincount)
    except Exception as e:
        print(e) 
        alpha = np.nan
    return alpha

meta = load_metadata_emerson(cohort=1, filtered=True, set_index=False)
meta['alpha'] = meta.apply(lambda s: alpha(emerson_processeddata_directory+'%s.tsv.gz'%s['Subject ID']), axis=1)
meta.to_csv('data/alpha_emerson_mincount%g.csv'%mincount)
