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
    counts = pd.read_csv(path, sep='\t')['count']
    alpha = mle_alpha_discrete(counts, cmin=mincount)
    return alpha

meta = load_metadata_britanova()
meta['alpha'] = meta.apply(lambda s: alpha(britanova_rawdata_directory+s['file_name']+'.gz'), axis=1)
meta.to_csv('data/alpha_britanova_mincount%g.csv'%mincount)
