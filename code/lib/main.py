import string
import numbers
import itertools
import sys, os
import numpy as np

np.seterr(all='raise')

def subsample(counts, n):
    """Randomly subsample from a vector of counts without replacement.

    Parameters
    ----------
    counts : Vector of counts (integers) to randomly subsample from.
    n : Number of items to subsample from `counts`. Must be less than or equal
        to the sum of `counts`.

    Returns
    -------
    Subsampled vector of counts where the sum of the elements equals `n`
    """
    n = int(n)
    unpacked = np.concatenate([np.repeat(np.array(i,), count) for i, count in enumerate(counts)])
    sample = np.random.choice(unpacked, size=n, replace=False)
    unique, counts = np.unique(sample, return_counts=True)
    return counts


def parametercheck(datadir, argv, paramscomb, nbatch):
    if not os.path.exists(datadir):
        print('datadir missing!')
        return False
    if not len(argv) > 1:
        if len(paramscomb) % nbatch != 0.0:
            print('incompatible nbatch', len(paramscomb), nbatch)
        else:
            print('number of jobs', len(paramscomb) // nbatch)
        return False
    return True

def params_combination(params):
    """Make a list of all combinations of the parameters."""
    # for itertools.product to work float entries have to be converted to 1-element lists
    params = [[p] if isinstance(p, numbers.Number) or isinstance(p, str) or hasattr(p, '__call__') else p for p in params]
    return list(itertools.product(*params))

def runcombinations_scalar(run, paramlist, paramnames, datadir='../data/', outname='out', nbatch=1, disp=True):
    """ Run the function `run` with all combinations of parameters"""
    paramscomb = params_combination(paramlist)
    if parametercheck(datadir, sys.argv, paramscomb, nbatch):
        njob = int(sys.argv[1])
        data = []
        for i in range(nbatch):
            n = (njob-1) * nbatch + i
            if disp:
                print(zip(paramnames[:len(paramscomb[n])], paramscomb[n]))
            res = run(*paramscomb[n])
            row = list(paramscomb[n])
            row.extend(res)
            data.append(row)
        np.savez_compressed(datadir + '%s%g' % (outname, njob), data=data, columns=paramnames)

def runcombinations_array(run, paramlist, paramnames, datadir='../data/',
    outname='out', disp=True, **kwargs):
    """ Run the function `run` which returns an array of variable length with all combinations of parameters"""
    paramscomb = params_combination(paramlist)
    print(paramscomb)
    print(sys.argv)
    if parametercheck(datadir, sys.argv, paramscomb, 1):
        n = int(sys.argv[1])
        if disp:
            print(list(zip(paramnames, paramscomb[n-1])))
        res = run(*paramscomb[n-1], **kwargs)
        np.savez_compressed(datadir + '%s%g' % (outname, n), res=res, paramnames=paramnames, params=paramscomb[n-1])
