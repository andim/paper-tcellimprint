import pandas as pd
import numpy as np
from lib import *

meta = load_metadata_emerson(filtered=True)

count_dict = {}
for SubjectID in meta.index:
    counts = pd.read_csv(emerson_processeddata_directory + SubjectID +'.tsv.gz')['counts']
    count_dict[SubjectID] = np.array(counts)

meta['totcount'] = meta.apply(lambda s: count_dict[s.name].sum(), axis=1)
meta['singletons'] = meta.apply(lambda s: np.sum(count_dict[s.name]==1), axis=1)
meta['largest'] = meta.apply(lambda s: count_dict[s.name][0], axis=1)
meta['10thlargest'] = meta.apply(lambda s: count_dict[s.name][10], axis=1)
meta['100thlargest'] = meta.apply(lambda s: count_dict[s.name][100], axis=1)
meta['1000thlargest'] = meta.apply(lambda s: count_dict[s.name][1000], axis=1)

meta.to_csv(data_directory + '/emerson-counts.csv', index=True)
