import pandas as pd
import numpy as np
from lib.config import *
import lib

meta = pd.read_csv(data_directory + 'metadata-britanova.txt', sep='\t')
mask = meta['age']>0
print(meta[mask].shape[0], meta.shape[0]-meta[mask].shape[0], meta[mask]['age'].min(), meta['age'].max())

Nsample = 5e5

count_dict = {}
for file_name in meta['file_name']:
    print(file_name)
    counts = pd.read_csv(britanova_rawdata_directory + file_name +'.gz', sep='\t', usecols=(0,))['count']
    counts = np.sort(counts) 
    count_dict[file_name] = counts

meta['totcount'] = meta.apply(lambda s: count_dict[s['file_name']].sum(), axis=1)
meta['singletons'] = meta.apply(lambda s: np.sum(count_dict[s['file_name']]==1), axis=1)
meta['largest'] = meta.apply(lambda s: count_dict[s['file_name']][-1], axis=1)
meta['10thlargest'] = meta.apply(lambda s: count_dict[s['file_name']][-10], axis=1)
meta['100thlargest'] = meta.apply(lambda s: count_dict[s['file_name']][-100], axis=1)
meta['1000thlargest'] = meta.apply(lambda s: count_dict[s['file_name']][-1000], axis=1)
meta['10000thlargest'] = meta.apply(lambda s: count_dict[s['file_name']][-10000], axis=1)

meta['sub_singletons'] = meta.apply(lambda s: np.sum(lib.subsample(count_dict[s['file_name']], Nsample)==1) if s['totcount'] > Nsample else np.nan, axis=1)

meta.to_csv(data_directory + '/britanova-counts.csv', index=False)
