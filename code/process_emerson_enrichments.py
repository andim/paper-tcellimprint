import pandas as pd
import numpy as np
from lib import *

metadata = load_metadata_adaptive_all(filtered=True)
print(metadata.shape)

zeroInsertions = {}
zeroInsertions_out = {}
out = {}
oneInsertions = {}
twoInsertions = {}
invdjdb = {}
invdjdb_lev1 = {}
bins = np.array([1, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000])
missing = []
for subjectid, meta in metadata.iterrows():
    print(subjectid)
    try:
        df = pd.read_csv(emerson_processeddata_directory + subjectid + '.tsv.gz') 
    except (FileNotFoundError, EOFError):
        missing.append(subjectid)
        print('missing')
        continue
    df['invdjdb_lev1'].fillna(False, inplace=True)
    df['oneInsertion'] = (df['totalInsertion'] == 1)
    df['twoInsertion'] = (df['totalInsertion'] == 2)
    df['zeroInsertion_out'] = df['zeroInsertion'] & (df['sequenceStatus'] != 'In')
    df['out'] = (df['sequenceStatus'] != 'In')
    dfg = df.groupby(pd.cut(df['rank'], bins-0.1))
    zeroInsertions[subjectid] = dfg.zeroInsertion.mean()
    oneInsertions[subjectid] = dfg.oneInsertion.mean()
    twoInsertions[subjectid] = dfg.twoInsertion.mean()
    invdjdb[subjectid] = dfg.invdjdb.mean()
    invdjdb_lev1[subjectid] = dfg.invdjdb_lev1.mean()
    zeroInsertions_out[subjectid] = dfg.zeroInsertion_out.mean()
    out[subjectid] = dfg.out.mean()

print('missing files', missing)
metadata = metadata.loc[~metadata.index.isin(missing)]

for i, rank in enumerate(bins[1:]):
    metadata['zeroInsertion%s'%rank] = [zeroInsertions[subjectid].iloc[i]
                                        for subjectid, row in metadata.iterrows()]
    metadata['zeroInsertion_out%s'%rank] = [zeroInsertions_out[subjectid].iloc[i]
                                        for subjectid, row in metadata.iterrows()]
    metadata['out%s'%rank] = [out[subjectid].iloc[i]
                                        for subjectid, row in metadata.iterrows()]
    metadata['oneInsertion%s'%rank] = [oneInsertions[subjectid].iloc[i]
                                        for subjectid, row in metadata.iterrows()]
    metadata['twoInsertion%s'%rank] = [twoInsertions[subjectid].iloc[i]
                                        for subjectid, row in metadata.iterrows()]
    metadata['invdjdb%s'%rank] = [invdjdb[subjectid].iloc[i]
                                        for subjectid, row in metadata.iterrows()]
    metadata['invdjdb_lev1_%s'%rank] = [invdjdb_lev1[subjectid].iloc[i]
                                        for subjectid, row in metadata.iterrows()]
metadata.to_csv(data_directory + '/emerson-enrichments.csv')
