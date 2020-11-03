import os.path
import glob
import numpy as np
import pandas as pd

from lib import *
from lib.pgen import calc_nt_pgen

files = sorted(glob.glob(emerson_rawdata_directory+'*.tsv.gz'))

vdjdb = pd.read_csv(data_directory +'vdjdb.slim.txt', sep='\t')
vdjdb = vdjdb[vdjdb['species'] == 'HomoSapiens']
vdjdb_aas = set(vdjdb['cdr3'])

if len(sys.argv) == 1:
    print(len(files))
else:
    filename = files[int(sys.argv[1])]
    outname = data_directory+'processed/emerson/' + filename.split('/')[-1]
    if not os.path.exists(outname):
        print(filename)
        df = pd.read_csv(filename, sep='\t')
        df.rename(columns=rename_adaptive1, inplace=True)
        # randomize order such that order within rank is random
        df = df.sample(frac=1).reset_index(drop=True)
        df['totalInsertion'] = df['n1Insertion']+df['n2Insertion'] 
        df['zeroInsertion'] = (df['n1Insertion'] == 0) & (df['n2Insertion'] == 0)
        df['rank'] = df['counts'].rank(ascending=False, method='first')
        df['invdjdb'] = df.apply(lambda row: (row['aa'] in vdjdb_aas), axis=1)
        df['invdjdb_lev1'] = df.apply(lambda row: dist1_levenshtein(row['aa'], vdjdb_aas), axis=1)
        df['pgen'] = df.apply(lambda row: calc_nt_pgen(row), axis=1)
        # order by rank
        df = df.sort_values('rank').reset_index(drop=True)
        df.to_csv(outname, index=False, columns=['rank', 'nt', 'aa', 'counts', 'sequenceStatus',
                                                 'zeroInsertion', 'totalInsertion',
                                                 'invdjdb', 'invdjdb_lev1', 'pgen'])
