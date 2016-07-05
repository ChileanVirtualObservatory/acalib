#!python

import numpy as np
cimport numpy as cnp
from clumpfind cimport AstKeyMap, cupidClumpFind


ctypedef cnp.int_t int_t
ctypedef cnp.float64_t float64_t
ctypedef cnp.ndarray ndarray

cdef float64_t[::1] _test():
	#clumpf( int type, int ndim, int *slbnd, int *subnd, 
	#void *ipd, double *ipv, double rms, AstKeyMap *config, int velax, 
	#int perspectrum, double beamcorr[ 3 ], int *backoff, int *status )
	cdef:
		int *slbnd
		int *subnd
		double *beamcorr
		int backoff = 0
		int status = 1
		ndarray[float64_t, ndim=3, mode='fortran'] arr
		AstKeyMap kmap 
	#slbnd[:] = [0,0,0]
	#subnd[:] = [10,10,10]
	#beamcorr = [0.,0.,0.]
	arr = np.random.random((10,10,10)).flatten(order='F')
	beamcorr = <double*> arr.data
	arr2 = cupidClumpFind(1, 3, slbnd, subnd, <double*> arr.data, NULL, 1., &kmap, 0, 0, beamcorr, &backoff, &status)
	return arr

def sum(a,b):
	return a+b

def test():
	return _test()