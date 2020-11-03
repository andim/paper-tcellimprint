from cycler import cycler

import sys
sys.path.append('..')
from lib import *

import re
import palettable

datadir = 'data/'

Nsample = 1e6

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

files = sorted(glob.glob(datadir + "earlydeath*.npz"), key=numericalSort)

def load(f):
    data = np.load(f)
    theta = float(data['params'][np.argwhere(data['paramnames']=='intro')[0]])
    C0 = float(data['params'][np.argwhere(data['paramnames']=='introsize')[0]])
    b0 = float(data['params'][np.argwhere(data['paramnames']=='birth')[0]])
    gamma = theta*C0/b0
    counts = data['res'][0]
    counts = np.copy(counts[counts>0])
    return gamma, counts

for i, f in enumerate(files):
    print(f)
    gamma, counts = load(f)
    counts = subsample(counts, int(Nsample)) 
    pd.DataFrame(counts).to_csv(datadir+'earlydeath_processed_gamma%g.csv.gz'%gamma, index=False)
