import numpy as np

aminoacids = 'ACDEFGHIKLMNPQRSTVWY'

def dist1_hamming(x, reference):
    """ Is the kmer x a Hamming distance 1 away from any of the kmers in the reference set"""
    for i in range(len(x)):
        for aa in aminoacids:
            if aa == x[i]:
                continue
            if x[:i]+aa+x[i+1:] in reference:
                return True
    return False

def dist1_levenshtein(x, reference):
    """ Is the kmer x a Levenshtein distance 1 away from any of the kmers in the reference set"""
    if not type(x) == type(''):
        return np.nan
    for i in range(len(x)):
        # deletion
        if x[:i]+x[i+1:] in reference:
            return True
        for aa in aminoacids:
            if aa == x[i]:
                continue
            # replacement
            if x[:i]+aa+x[i+1:] in reference:
                return True
            # insertion
            if x[:i]+aa+x[i:] in reference:
                return True
    return False
