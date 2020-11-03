# library for explicit antigen simulation
import numpy as np


def get_rates(c, a, K, t, nu=lambda t: 1.0, mu=lambda t: 1.0, b=lambda S: 1.0,
              d=lambda S: 1.0, availability_function=None):
    """ Implements the following dynamical equations
    dci/dt = (b(Si) * nu(t) - d(Si) * mu(t))*ci
    with S_i = K_ij a_j g(T_j) and T = c K

    T: antigenic trigger
    """
    # compute stimulus: interaction mat. x ag. conc. x competition
    if availability_function is not None:
        # antigenic trigger vector (1 entry per clone)
        antigenic_trigger = c.dot(K)
        S = K.dot(a * availability_function(antigenic_trigger))
    # no competition function
    else:
        # stimulus for clone i= interaction matrix
        S = K.dot(a)

    # returns absolute rate, total rate
    bS, dS = b(S) * nu(t), d(S) * mu(t)
    return bS, dS


def sample_discrete(probs):
    """
    Randomly sample an index with probability given by probs.
    """
    q = np.random.rand()
    i = 0
    p_sum = 0.0
    while p_sum < q:
        p_sum += probs[i]
        i += 1
    return i - 1

def deterministic_evolution_steps(c, a, K, t, dt, theta_c, c_cutoff,
                                  rate_function, prng=np.random, stochasticity='no'):
    """
    Draw time until next introduction event from Poisson process, based on
    introduction rates. Loop as long as no new introduction occurs.
    Stochasticity:  'no': fully deterministic,
                    'linearized': deterministic rates modified by birth-death
                        noise
    """
    # determine next intro event
    t_old = t
    t_delta = prng.exponential(1.0 / theta_c)

    while t < t_old + t_delta:
        # make time step fit next intro event time
        if t + dt > t_old + t_delta:
            thisdt = t_old + t_delta - t
        else:
            thisdt = dt

        # adapt rates
        birth_c, death_c = rate_function(c, a, K, t)

        # evolve clone sizes in the continuum approximation (deterministic)
        if stochasticity == 'linearized':  # with birth-death-noise
            c += (birth_c - death_c) * c * thisdt + ((birth_c + death_c) * c * thisdt) ** .5 * prng.normal()
        else:
            c += (birth_c - death_c) * c * thisdt
        t += thisdt

        # clones / antigens lower than threshold are zeroed
        c[c < c_cutoff] = 0.0

    return c, t_old + t_delta


def discrete_stochastic_evolution_step(c, a, K, t, dt, theta_c,
                                       rate_function, prng=np.random, verbose=False):
    """Simulate clonal evolution explicitly using discrete probabilistic events,
    i.e. the immigration (introduction), birth and death of clones.
    """
    # determine rates
    birth_c, death_c = rate_function(c, a, K, t)
    birth_c = birth_c * c
    death_c = death_c * c
    total_rate = theta_c + np.sum(birth_c) + np.sum(death_c)

    # draw time until next event
    dt = prng.exponential(1.0 / total_rate)
    if verbose:
        print("dt=\t%f" % dt)
    t += dt

    # assemble probability space
    probs = np.zeros(2 * len(birth_c) + 1)
    probs[0] = theta_c / total_rate
    probs[1:len(birth_c) + 1] = birth_c / total_rate
    probs[1 + len(birth_c):] = death_c / total_rate
    if verbose:
        print('Probability space:', probs)

    # decide event
    ind = sample_discrete(probs)
    intro_bool = False
    if ind == 0:
        intro_bool = True
        if verbose:
            print('Introduction event')
    elif ind <= len(birth_c):
        c[ind - 1] += 1
        if verbose:
            print('Proliferation event')
    else:
        c[ind - 1 - len(birth_c)] -= 1
        if verbose:
            print('Death event')
    return c, t, intro_bool


def execute_introduction_event(c, a, c_new, K, p_binding, c_cutoff,
                               prng=np.random, verbose=False):
    """Adapt the clonal array c and the affinity matrix K when a new clone
    enters the repertoire.
    """
    indmin = np.argmin(c)
    replace = c[indmin] < c_cutoff
    # replace existing (too) small clone or add to rep.
    if replace:
        c[indmin] = c_new
        if verbose:
            print('Replacement! Clone nr. %g' % indmin)
    else:
        c = np.concatenate((c, (c_new,)))
        if verbose:
            print('Addition! Clone nr. %s' % str(len(c) - 1))

    # handle affinity matrix
    newaffinity = prng.random_sample(len(a)) < p_binding
    if replace:  # replace respective existing binding prob.
        K[indmin, :] = newaffinity
    else:  # or, if new clone, add column to K with new entries
        K = np.concatenate((K, newaffinity[None, :]), axis=0)
    return c, K


def print_progress(c, t, K, t_end, flag):
    """Utility to print simulation progress with some key
    debugging information.
    """
    if (t > (t_end / 20.0)) and flag == True:
        print('1/20 of the simulation is done.')
        flag = False

    # debugging information
    print('\n======================\nTime: %g' % (t))
    print('Clones:\n', c)
    print('\nInteractions in Matrix K:\n', np.array(K))
    print('Dimension of Matrix K:\t%s' % str(K.shape))
    return flag


def simulate(rate_function=None, theta_c=1.0,
             t_start=0, t_end=1, dt=0.01,
             a_initial=0.0, c_initial=0.0,
             Nc=1, Na=1,
             a_new=1.0, c_new=1.0,
             p_binding=None,
             c_cutoff=0.1,
             stochasticity='no',
             verbose=False, prng=np.random, seed=1996,
             ):
    """Simulates explicit competition of T cells for antigens in a discrete,
    deterministic or linearized (deterministic with birth-death noise) way.

    rate_func: function taking c, a, K and returning rates of change for c, a
    theta_c: introduction (immigration) rate of new clones
    t_start, t_end: beginning and end of simulation
    dt: simulation time step
    a_initial, c_initial: antigen / clone size at beginning of simulation
    Na, Nc: initial number of different antigens / clones
    a_new, c_new: antigen / clone size of newly introduced clones during the sim
    p_binding: probability for a randomly drawn clone to bind to a given antigen
        The affinity matrix K is constructed by drawing random bindings
        according to this prob.
    c_cutoff: threshold clone size below which a clone is considered as extinct
        and gets removed or replaced
    stochasticity: one of 'no', 'exact', 'linearized'
    verbose: enable printing of debugging information
    prng: pseudo random number generator
    seed: seed to initialize prng
    """
    if verbose:
        print(
            "### INIT\ntheta_c=\t%g\na_new=\t\t%g\nc_new=\t\t%g\nNa=\t\t%g\nNc=\t\t%g\np_binding=\t%g\nt_end=\t\t%g\nc_cutoff=\t%g\n" % (
                theta_c, a_new, c_new, Na, Nc, p_binding, t_end, c_cutoff))

    # initialize arrays
    prng.seed(seed=seed)
    c = np.ones(Nc) * c_initial  # fill initial a/c arrays
    a = np.ones(Na) * a_initial
    K = prng.random_sample((Nc, Na)) < p_binding  # generate interaction matrix

    # main loop
    t = t_start
    progressflag = True

    while t < t_end:
        introevent = False

        # measure progress
        if verbose:
            progressflag = print_progress(c, t, K, t_end, progressflag)

        # simulate until introduction event occurs
        if stochasticity in ['no', 'linearized']:
            c, t = deterministic_evolution_steps(c, a, K, t, dt, theta_c, c_cutoff,
                                                 rate_function, prng=prng, stochasticity=stochasticity)
            introevent = True

        # simulate one discrete time step with one event
        elif stochasticity is 'discrete':
            c, t, introevent = discrete_stochastic_evolution_step(c, a, K, t, dt,
                                                                  theta_c, rate_function, prng=prng, verbose=verbose)
        else:
            raise Exception('Unknown stochasticity!')

        if introevent:
            c, K = execute_introduction_event(c, a, c_new, K, p_binding, c_cutoff,
                                              prng=prng, verbose=verbose)

    return c, a, K
