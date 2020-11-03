import sys
sys.path.append('..')
import pandas as pd
import numpy as np
from lib import *

metadata = pd.read_csv('data/meta.csv', index_col=0)
print(metadata)

early_probability = 0.07
background_probability = 0.02
early_age = 0.05
bins = np.array([1, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000])#, 200000, 500000, 1000000])

zeroInsertions = {}
missing = []
for subjectid, meta in metadata.iterrows():
    print(subjectid)
    try:
        df = pd.read_csv('data/ff_' + str(subjectid) + '.csv.gz') 
    except (FileNotFoundError, EOFError):
        missing.append(subjectid)
        print('missing')
        continue
    #df['early'] = df['Age'] > (meta['Age']-early_age)
    #df['zeroInsertion'] = np.random.binomial(1, df['early']*max_probability + background_probability)
    ts = df['Age']
    fraction_early = early_age/(ts+1e-10) # add epsilon to avoid divide by zero
    zeroProb = np.where(ts<early_age, early_probability,
                     background_probability*(1-fraction_early) + early_probability*fraction_early)
    df['zeroInsertion'] = np.random.binomial(1, zeroProb)
    print(df['zeroInsertion'].mean(), meta['Age'])
    # randomize order such that order within rank is random
    df = df.sample(frac=1).reset_index(drop=True)
    df['rank'] = df['counts'].rank(ascending=False, method='first')
 
    dfg = df.groupby(pd.cut(df['rank'], bins-0.1))
    zeroInsertions[subjectid] = dfg.zeroInsertion.mean()

print('missing files', missing)
metadata = metadata.loc[~metadata.index.isin(missing)]

for i, rank in enumerate(bins[1:]):
    metadata['zeroInsertion%s'%rank] = [zeroInsertions[subjectid].iloc[i]
                                        for subjectid, row in metadata.iterrows()]
metadata.to_csv('data/enrichments.csv')
