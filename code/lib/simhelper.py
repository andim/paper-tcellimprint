from . import *
import lib.simexplag as sim

debug = False # debugging mode

def run(birth_ini, death_ini, intro, introsize, tlim, dyn, sigma=0.0,
    debug=False, **kwargs):
    """tlim: maximum time for simulation and clone size recording,
    sigma: flucuating fitness standard deviation
    introgen: generative function for drawing introduction sizes.
        (default: constant value specified in 'introsize')
    """
    prng = np.random
    if "seed" in kwargs:
        prng.seed(seed=kwargs.get("seed"))

    if "introgen" not in kwargs:  # default: constant introduction size
        introgen = lambda: introsize
    else:
        introgen = kwargs.get("introgen")
        
    print("\nb=%f, d=%f, i=%f, C0=%f, dyn=%s" %
        (birth_ini, death_ini, intro, introsize, dyn))

    def birthf(t):
        if death_ini == 0:
            return birth_ini / ((birth_ini + intro*introsize)*t)
        return death_ini / ((1.0 + intro*introsize/birth_ini) * (1.0-np.exp(-death_ini * t)))
    def deathf(t):
        return death_ini

    cumulrun = intro * tlim
    cr   = 0
    clonesizes = []
    cloneages = []
    while cr < cumulrun:
        cr += 1
        tstart = prng.random() * tlim
        t = tstart
        cloneages.append(tlim - tstart)
        clonesize = introgen()

        while True:

            if debug:
                print("Time:\t\t%f / %f" % (t, tlim))
                print("Cumulrun:\t%d / %d" % (cr+1, cumulrun))
                print("Adapt. deathr.:\t%f" % deathf(t))
                print("Adapt. birthr.:\t%f" % birthf(t))


            if sigma and ((birthf(t) + deathf(t))*5 < clonesize):#((birthf(t) + deathf(t))/clonesize < sigma**2):
                f0 = -sigma**2/2 + birthf(t) - deathf(t)
                deltat = 0.01/birthf(t) #np.abs(f0)/10.0
                finished = False
                if t + deltat > tlim:
                    deltat = tlim-t
                    finished = True
                factor = prng.lognormal(mean=deltat*f0,
                                sigma=sigma*np.sqrt(deltat))
                clonesize = int(round(clonesize*factor))
                t += deltat
                if finished:
                    break
            else:
                rr_tot = lambda t: (birthf(t) + deathf(t))*clonesize
                t = next_event(t, rr_tot, dyn, tlim, debug)
                if t > tlim:
                    break
                r = prng.random()
                if r < (birthf(t)*clonesize)/rr_tot(t):
                    clonesize += 1
                else:
                    clonesize -= 1
            if clonesize == 0:
                break
        clonesizes.append(clonesize)
        #if cr % (cumulrun//10) == 0:
        #    print('%g of %g'%(cr, cumulrun))
    if debug:
        print("\n+++RESULTS+++\n")
        print(clonesizes, cloneages)

    return np.array([clonesizes, cloneages])

def run_exact_const(birth_ini, death_ini, intro_ini, *args, **kwargs):
    print("\nb=%f, d=%f, i=%f" % (birth_ini, death_ini, intro_ini))
    def birth(c, t):
        return birth_ini*c/np.sum(c)
    def death(c, t):
        return death_ini*c
    def intro(c, t):
        return intro_ini
    return run_exact(birth, death, intro, *args, **kwargs)

def run_exact_K(birth_ini, death_ini, intro_ini, K, *args, **kwargs):
    print("\nb=%f, d=%f, i=%f, K=%f" % (birth_ini, death_ini, intro_ini, K))
    def birth(c, t):
        return birth_ini*c/(np.sum(c) + K)
    def death(c, t):
        return death_ini*c
    def intro(c, t):
        return intro_ini
    return run_exact(birth, death, intro, *args, **kwargs)

def run_exact(birth, death, intro, introsize, tlim, **kwargs):
    "tlim: maximum time for simulation and clone size recording"
    prng = np.random
    if "seed" in kwargs:
        prng.seed(seed=kwargs.get("seed"))
    if "introgen" not in kwargs:  # default: constant introduction size
        introgen = lambda: introsize
    else:
        introgen = kwargs.get("introgen")

    t = prng.exponential(1.0/intro(np.array([]), 0.0))
    if debug:
        print("t=\t%e" % t)
    c = np.array([1])

    # initialize clone size recording
    if "recordingnumber" in kwargs:
        recordingnumber = kwargs.get("recordingnumber")
    else:
        recordingnumber = None
    if recordingnumber != None:
        clonesizes_timeseries = sim.init_timeseries_recording(recordingnumber=recordingnumber)

    while t < tlim:
        birth_c = birth(c, t)
        death_c = death(c, t)
        intro_c = intro(c, t)
        totrate = intro_c + np.sum(birth_c) + np.sum(death_c)
        dt = prng.exponential(1.0/totrate)
        if debug:
            print("dt=\t%e" % dt)
        t += dt
        probs = np.zeros(2*len(birth_c)+1)
        probs[0] = intro_c/totrate
        probs[1:len(birth_c)+1] = birth_c/totrate
        probs[1+len(birth_c):] = death_c/totrate
        if debug:
            print('Probability space:', probs)
        ind = sim.sample_discrete(probs)
        if ind == 0:
            if debug:
                print('Introduction event')
            indmin = np.argmin(c)
            if c[indmin] < 1:
                if "introdistr" not in kwargs:
                    c[indmin] = introsize
                else:
                    c[indmin] = introgen()
                if debug:
                    print('Replacement! Clone nr. %g' % indmin)
            else:
                c = np.concatenate((c, (introsize,)))
                if debug:
                    print('Addition! Clone nr. %s' % str(len(c)))
        elif ind <= len(birth_c):
           c[ind-1] += 1
           if debug:
               print('Proliferation event')
        else:
            c[ind-1-len(birth_c)] -= 1
            if debug:
                print('Death event')

        # clone sizes recording
        if (recordingnumber != None and c[0] != 0.0):
        # do not record dummy clone with size 0 in the beginning
            clonesizes_timeseries=sim.record_timeseries(clonesizes_timeseries,
                recordingnumber, c, t)

        # exit condition
        if recordingnumber != None and len(c) >= recordingnumber:
            break

    if recordingnumber == None:
        return np.array([c, np.zeros_like(c)])
    else:
        return np.array([c, clonesizes_timeseries])
