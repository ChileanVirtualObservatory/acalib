import cython

import numpy as np
cimport numpy as np

cdef extern:
    void diff_c(double* cumPixels, double* diff, int n)
    void seg_c (double* diff, double* boxing, int n)
    void eros_c (double* boxing, double* blocking, int n)

@cython.boundscheck(False)
@cython.wraparound(False)
def diff(np.ndarray[double,ndim=1, mode="c"] cumPixels):
    cdef np.int32_t n = len(cumPixels)
    cdef np.ndarray[double, ndim=1] diff = np.zeros(n,dtype=np.float64)
    diff_c(<double *> &cumPixels[0], <double *> &diff[0],n)
    
    return diff

@cython.boundscheck(False)
@cython.wraparound(False)
def seg(np.ndarray[double,ndim=1, mode="c"] diff):
    cdef np.int32_t n = len(diff)
    cdef np.ndarray[double, ndim=1] boxing = np.zeros(n,dtype=np.float64)
    seg_c(<double *> &diff[0], <double *> &boxing[0], n)

    return boxing

@cython.boundscheck(False)
@cython.wraparound(False)
def eros(np.ndarray[double,ndim=1, mode="c"] boxing):
    cdef np.int32_t n = len(boxing)
    cdef np.ndarray[double, ndim=1] blocking = np.zeros(n,dtype=np.float64)
    eros_c(<double *> &boxing[0], <double *> &blocking[0], n)

    return boxing    
