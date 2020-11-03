import sys
sys.path.append('..')
from lib import *
from lib.simhelper import *

np.seterr(all='raise')

# initial population  and time parameters
intros    = [2e3] # introduction / innovation rate
introsizes  = [1] # clone size at introduction
tlims = [5.0] # end time
births    = [2e4] # birth rate
deaths    = [0.2] # death rate
Ks        = [0.0, 10.0, 100.0, 1000.0]

runcombinations_array(run_exact_K, (births, deaths, intros, Ks, introsizes, tlims),
        ('birth', 'death', 'intro', 'K', 'introsize', 'tlim'), datadir='data/', outname='K')
