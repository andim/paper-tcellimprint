import os.path
import glob
import numpy as np
import pandas as pd

from lib import *
from lib.pgen import calc_nt_pgen

files = pd.read_csv(data_directory + 'metadata-lindau.csv', index_col=0)
print(files)

vdjdb = pd.read_csv(data_directory +'vdjdb.slim.txt', sep='\t')
vdjdb = vdjdb[vdjdb['species'] == 'HomoSapiens']
vdjdb_aas = set(vdjdb['cdr3'])

for subjectid in files.index:
    print(subjectid)
    outname = data_directory+'processed/emerson/' + subjectid + '.tsv.gz' 
    if not os.path.exists(outname):
        filename = data_directory + 'raw/lindau/' + subjectid[-2:] + '_PBMC.tsv.gz'
        df = pd.read_csv(filename, sep='\t')
        df.rename(columns=rename_adaptive2, inplace=True)
        # randomize order such that order within rank is random
        df = df.sample(frac=1).reset_index(drop=True)
        df['totalInsertion'] = df['n1Insertion']+df['n2Insertion'] 
        df['zeroInsertion'] = (df['n1Insertion'] == 0) & (df['n2Insertion'] == 0)
        df['rank'] = df['counts'].rank(ascending=False, method='first')
        df['invdjdb'] = df.apply(lambda row: (row['aa'] in vdjdb_aas), axis=1)
        df['invdjdb_lev1'] = df.apply(lambda row: dist1_levenshtein(row['aa'], vdjdb_aas), axis=1)
        # different nucleotide format
        # does not allow the same trimming used elsewhere to find CDR3 sequence
        # therefore set all pgens to nan
        df['pgen'] = df.apply(lambda row: np.nan)
        # order by rank
        df = df.sort_values('rank').reset_index(drop=True)
        df.to_csv(outname, index=False, columns=['rank', 'nt', 'aa', 'counts', 'sequenceStatus',
                                                 'zeroInsertion', 'totalInsertion',
                                                 'invdjdb', 'invdjdb_lev1', 'pgen'])
