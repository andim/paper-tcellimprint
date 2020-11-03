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
dyn       = ["poisson_thinned"] # clone dynamics

runcombinations_array(run, (births, deaths, intros, introsizes, tlims, dyn ), ('birth', 'death', 'intro', 'introsize', 'tlim', 'dyn'), datadir='data/', outname='mf')
