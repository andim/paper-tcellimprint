import pandas as pd
import numpy as np
from lib import *

meta = load_metadata_emerson(cohort=1, filtered=False)
print(meta)

print('total number of datasets', meta.shape[0])
print(list(meta[meta['Age'].isna()].index))
meta = meta[~meta['Age'].isna()]
print('filter out datasets with missing age', meta.shape[0])

count_dict = {}
mask = []
for SubjectID in meta.index:
    print('.', end='', flush=True)
    df = pd.read_csv(emerson_rawdata_directory + SubjectID +'.tsv.gz', nrows=1, sep='\t')
    mask.append(df['templates'].isnull())
mask = np.array(mask)
print(list(meta[mask].index))
meta = meta[~mask]
print('filter out datasets with missing template assignment', meta.shape[0])
meta.to_csv(data_directory + 'metadata-emerson-filtered.csv')

