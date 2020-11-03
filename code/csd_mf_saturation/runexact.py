import sys
sys.path.append('..')
from lib import *
from lib.simhelper import *

np.seterr(all='raise')

# initial population  and time parameters
intros    = [2e3] # introduction / innovation rate
introsizes  = [1, 1, 1] # clone size at introduction
tlims = [5.0] # end time
births    = [2e4] # birth rate
deaths    = [0.2] # death rate


runcombinations_array(run_exact_const, (births, deaths, intros, introsizes, tlims),
        ('birth', 'death', 'intro', 'introsize', 'tlim'), datadir='data/', outname='exact')
