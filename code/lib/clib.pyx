import numpy as np
cimport numpy as np

def sample_discrete(np.ndarray[np.double_t, ndim=1] probs):
    """
    Randomly sample an index with probability given by probs.
    """
    cdef double q = np.random.rand()
    cdef int i = 0
    cdef double p_sum = 0.0
    while p_sum < q:
        p_sum += probs[i]
        i += 1
    return i - 1


