import pandas as pd
import numpy as np
from lib import *

metadata = load_metadata_adaptive_all(filtered=True)
print(metadata.shape)

zeroInsertions = {}
zeroInsertionsCells = {}
missing = []
for subjectid, meta in metadata.iterrows():
    print(subjectid)
    try:
        df = pd.read_csv(emerson_processeddata_directory + subjectid + '.tsv.gz', usecols=(3,5)) 
    except (FileNotFoundError, EOFError):
        missing.append(subjectid)
        print('missing')
        continue
    zeroInsertions[subjectid] = df.zeroInsertion.mean()
    weights = df['counts']
    weights /= np.sum(weights)
    zeroInsertionsCells[subjectid] = np.average(df.zeroInsertion, weights=weights)

print('missing files', missing)
metadata = metadata.loc[~metadata.index.isin(missing)]
metadata['zeroInsertions'] = [zeroInsertions[subjectid] for subjectid, row in metadata.iterrows()]
metadata['zeroInsertionsCells'] = [zeroInsertionsCells[subjectid] for subjectid, row in metadata.iterrows()]
metadata.to_csv(data_directory + '/emerson-totals.csv')
