from functools import reduce
import pandas as pd
import numpy as np

from .config import *

def multimerge(dfs, on, suffixes=None):
    "Merge multiple dataframes with suffix support"
    if suffixes:
        dfs_new = []
        for df, suffix in zip(dfs, suffixes):
            df = df.set_index(on)
            dfs_new.append(df.add_suffix('_'+suffix))
        return reduce(lambda left, right: pd.merge(left, right,
                                                   right_index=True, left_index=True,
                                                   how='outer'),
                      dfs_new)
    return reduce(lambda left, right: pd.merge(left, right, on), dfs)

metadata_columns = ['Age', 'Sex', 'Known CMV status', 'CMVpos']

metadata_kwargs = dict(na_values=('Unknown', 'ND'),
                       true_values=('+','True'),
                       false_values=('-','False'))

def load_metadata_britanova():
    meta = pd.read_csv(data_directory + 'metadata-britanova.txt', sep='\t')
    return meta

def load_metadata_adaptive_all(**kwargs):
    """ Load metadata from both Emerson cohorts and the Lindau mini-cohort. """
    kwargs['cohort'] = 'both'
    metadata_emerson = load_metadata_emerson(**kwargs)
    metadata_lindau = load_metadata_lindau(**kwargs)
    metadata = pd.concat([metadata_emerson, metadata_lindau])
    return metadata

def load_metadata_emerson(set_index=True, cohort='both', filtered=True):
    """"
    Load metadata of the Emerson Nat. Gen. study
    set_index = use Subject ID as index
    cohort = 'both', 1, 2
    """
    filename_emerson = 'metadata-emerson%s.csv' % ('-filtered' if filtered else '')
    columns = metadata_columns.copy()
    kwargs = metadata_kwargs.copy()
    if set_index:
        kwargs['index_col'] = 0
    else:
        columns.append('Subject ID')

    # filter out extra columns and rename CMV column
    def process(meta):
        meta = meta.filter(columns)
        return meta.rename(columns={'Known CMV status' : 'CMVpos'})

    if cohort == 1:
        meta = pd.read_csv(data_directory + filename_emerson, **kwargs)
        meta = process(meta)
    elif cohort== 2:
        meta = pd.read_csv(data_directory + 'metadata-emerson2.csv', **kwargs)
        meta = process(meta)
    elif cohort == 'both':
        meta1 = pd.read_csv(data_directory + filename_emerson, **kwargs)
        meta1 = process(meta1)
        meta2 = pd.read_csv(data_directory + 'metadata-emerson2.csv', **kwargs)
        meta2 = process(meta2)
        meta = pd.concat([meta1, meta2])
    else:
        raise ValueError(cohort)
    # fix some inconsistent naming
    if set_index:
        # there is no file HIP08203, but there is a file HIP08200 with matching metadata
        meta.rename(index=dict(HIP08203='HIP08200'), inplace=True)
        # the second cohort filenames end with _MC1
        meta.rename(index=lambda x: x + '_MC1' if (x.startswith('Keck') and not x.endswith('_MC1')) else x, inplace=True)
    return meta

def load_metadata_lindau(set_index=True, **kwargs):
    """"
    Load metadata of the Lindau J Immunol study
    set_index = use Subject ID as index
    """
    kwargs = metadata_kwargs.copy()
    columns = metadata_columns.copy()
    if set_index:
        kwargs['index_col'] = 0
    else:
        columns.append('Subject ID')
    meta = pd.read_csv(data_directory + 'metadata-lindau.csv', **kwargs)
    meta = meta.filter(columns)
    meta.rename(columns={'Known CMV status' : 'CMVpos'}, inplace=True)
    return meta

## column name standardization
# adaptive column names change between different files
# we thus provide a standardized mapping of column names

# used e.g. for the Emerson cohorts (when downloaded as standalone zip file)
rename_adaptive1 = dict(n1_insertions='n1Insertion',
              n2_insertions='n2Insertion',
              rearrangement='nt',
              amino_acid='aa',
              templates='counts',
              frame_type='sequenceStatus',
            )
 
# used e.g. in the Lindau study (when downloaded from immunoaccess)
rename_adaptive2 = dict(n1_insertions='n1Insertion',
              n2_insertions='n2Insertion',
              nucleotide='nt',
              aminoAcid='aa',
              sequenceStatus='sequenceStatus',
            )
rename_adaptive2['count (templates/reads)'] = 'counts'


def load_adaptive_data(filename, rank=False, freqs=False,
        unproductive=False, zeroInsertion=False):
    usecols = [0, 2]
    names = ['nucleotide', 'counts']
    if unproductive:
        usecols.append(1)
        names.append('aminoAcid')
    if zeroInsertion:
        usecols.extend([27, 30])
        names.extend(['n1Insertion', 'n2Insertion'])
    df = pd.read_csv(filename,
                     sep='\t', header=None,
                     usecols=usecols, names=names,
                     skiprows=1
                     )
    if freqs:
        df['freqs'] = df['counts']/df['counts'].sum()
    if rank:
        df['rank'] = df['counts'].rank(ascending=False, method='first')
    if unproductive:
        df['unproductive'] = ~df['aminoAcid'].notnull()
    if zeroInsertion:
        df['zeroInsertion'] = (df['n1Insertion'] == 0) & (df['n2Insertion'] == 0)
    return df

def load_chu(subject, phenotype, **kwargs):
    meta = pd.read_csv(data_directory+'metadata-chu.csv', parse_dates=['date'])
    selected_datasets = meta[(meta['subject'] == subject) & (meta['phenotype'] == phenotype)].sort_values('date')

    dfs = []
    for ind, row in selected_datasets.iterrows():
        df = load_adaptive_data(chu_rawdata_directory+row['path'], freqs=True, **kwargs)
        dfs.append(df)

    time = np.array((selected_datasets['date'].iloc[:]-selected_datasets['date'].iloc[0]).dt.days.values, dtype=np.float)
    time /= 365.

    dfm = multimerge(dfs, 'nucleotide', suffixes=[str(i) for i in range(len(dfs))])
    return time, dfm

