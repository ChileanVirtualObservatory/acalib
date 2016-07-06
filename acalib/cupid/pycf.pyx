#!python

import numpy as np
cimport numpy as cnp
from pycf cimport cupidClumpFind
#from clumpfind cimport AstKeyMap, cupidClumpFind


ctypedef cnp.int_t int_t
ctypedef cnp.float64_t float64_t
ctypedef cnp.ndarray ndarray

cdef int[:] _clumpfind(ndarray[double, ndim=1, mode="c"] data,config,rms,ndarray[int, ndim=1, mode="c"] shape):
	#clumpf( int type, int ndim, int *slbnd, int *subnd, 
	#void *ipd, double *ipv, double rms, AstKeyMap *config, int velax, 
	#int perspectrum, double beamcorr[ 3 ], int *backoff, int *status )
	cdef:
		int *slbnd = [0,0,0]
 		#ndarray[float64_t, ndim=1, mode='fortran'] inp_arr
		#int *s
		double *beamcorr = [0.,0.,0.]
		int backoff = 0
		int status = 0
#		ndarray[float64_t, ndim=1, mode='fortran'] inp_arr
		int [:] out_arr
	cdef PyObject* kmap = <PyObject *> config
	#inp_arr = np.random.random((10,10,10)).flatten(order='F')
	#inp_arr = data.flatten(order='F')
	elm = data.size
	out_arr = <int[:elm]> cupidClumpFind(1, 3,slbnd,&shape[0],&data[0], NULL, rms, kmap, 0, 0, beamcorr, &backoff, &status)
	#out_arr = .asarray(out_arr)
	return out_arr

#TODO Generalize to ndims!!
def clumpfind(data,config,rms):
	mv = _clumpfind(data.flatten(order='F'),config,rms,np.asarray(data.shape,dtype=np.int32))
	cb = np.reshape(mv,data.shape,order='F')
	return cb
