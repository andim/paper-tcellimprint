import sys
sys.path.append('..')
import lib.simexplag as simexplag
from lib import runcombinations_array

# fixed parameter set
mu_0 = 0.2 * 1e-3  # death rate
nu_0 = 10.0  # birth rate

a_new = 1.0  # antigen introduction size during sim
a_initial = 1.0  # antigen sizes at beginning of sim

theta_c = 1.0  # clonal introduction rate
c_new = 1.0  # clonal introduction size during simulate
c_initial = 1.0  # clone sizes at beginning of sim

# parameter combinations set
t_end = [10000.0]  # simulation time (units set such that theta_c = 1)
Na = [1000]  # number of antigens
p_binding = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]  # clone-antigen binding probability


def run(Na, p_binding, t_end):
    """Main run function, returns a clone size array at the end of the simulation.
    """
    # (in general) time dependent rates
    nu = lambda t: nu_0 / Na
    mu = lambda t: mu_0
    avail_func = lambda T: 1.0 / (1.0 + T)
    rate_func = lambda c, a, K, t: simexplag.get_rates(
        c, a, K, t, nu=nu, mu=mu,
        availability_function=avail_func,
        b=lambda S: S, d=lambda S: 1.0)
    c, a, K = simexplag.simulate(
        c_new=c_new, a_new=a_new, c_cutoff=0.0,
        dt=0.01, stochasticity='no',
        a_initial=a_initial, c_initial=c_initial,
        verbose=False, p_binding=p_binding,
        rate_function=rate_func, t_end=t_end,
        theta_c=theta_c, Na=Na)
    return c


runcombinations_array(run, (Na, p_binding, t_end), ('Na', 'pb', 'tend'), datadir='data/', outname='sim')
