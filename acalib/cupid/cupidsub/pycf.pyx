#!python

import numpy as np
cimport numpy as cnp
from clumpfind cimport cupidClumpFind

ctypedef cnp.float64_t float64_t

cdef float64_t[::1] _test():
	#clumpf( int type, int ndim, int *slbnd, int *subnd, 
	#void *ipd, double *ipv, double rms, AstKeyMap *config, int velax, 
	#int perspectrum, double beamcorr[ 3 ], int *backoff, int *status )
	cdef:
		int *slbnd = [0,0,0]
		int *subnd = [10,10,10]
		double *beamcorr = [0,0,0]
		int *backoff = [0]
		int *status = [1]
	arr = np.random.random((10,10,10)).flatten(order='F')
	arr2 = cupidClumpFind(1, 3, slbnd, subnd, <double*> arr.data, NULL, 1., dict(), 0, 0, beamcorr, backoff, status)
	return arr

def sum(a,b):
	return a+b

def test():
	return _test()