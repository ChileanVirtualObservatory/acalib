#!python
import numpy as np
cimport numpy as cnp
from pycupid cimport cupidClumpFind, cupidFellWalker

ctypedef cnp.ndarray ndarray


cdef int[:] _clumpfind(ndarray[double, ndim=1, mode="c"] data, config, 
						double rms, ndarray[int, ndim=1, mode="c"] shape):
	cdef:
		#Dependent parameters of ndim
		int ndim = shape.size
		ndarray[int, ndim=1, mode="c"] _slbnd = np.zeros(ndim, dtype=np.int32)
		ndarray[double, ndim=1, mode="c"] _beamcorr = np.zeros(ndim)
	
	cdef:
		#Parameters of cupidClumpFind
		int *slbnd = &_slbnd[0]
		int *subnd = &shape[0]
		void *ipd = &data[0]
		double *ipv = NULL
		PyObject* kmap = <PyObject *> config
		int velax = 0
		int perspectrum = 0
		double *beamcorr = &_beamcorr[0]
		int backoff = 0
		int status = 0

	cdef:
		#Output array
		int elm = data.size
		int[:] out_arr

	#Main function call
	out_arr = <int[:elm]> cupidClumpFind(1, ndim, slbnd, subnd, ipd, ipv, rms, kmap, velax, 
										perspectrum, beamcorr, &backoff, &status)
	return out_arr


cdef int[:] _fellwalker(ndarray[double, ndim=1, mode="c"] data, config, 
						double rms, ndarray[int, ndim=1, mode="c"] shape):
	cdef:
		#Dependent parameters of ndim
		int ndim = shape.size
		ndarray[int, ndim=1, mode="c"] _slbnd = np.zeros(ndim, dtype=np.int32)
		ndarray[double, ndim=1, mode="c"] _beamcorr = np.zeros(ndim)
	
	cdef:
		#Parameters of cupidClumpFind
		int *slbnd = &_slbnd[0]
		int *subnd = &shape[0]
		void *ipd = &data[0]
		double *ipv = NULL
		PyObject* kmap = <PyObject *> config
		int velax = 0
		int perspectrum = 0
		double *beamcorr = &_beamcorr[0]
		int status = 0

	cdef:
		#Output array
		int elm = data.size
		int[:] out_arr

	#Main function call
	out_arr = <int[:elm]> cupidFellWalker(1, ndim, slbnd, subnd, ipd, ipv, rms, kmap, velax, 
										perspectrum, beamcorr, &status)
	return out_arr



def clumpfind(data not None, config not None, rms):
	data = data.copy()
	shape = np.asarray(data.shape,dtype=np.int32)-1
	mv = _clumpfind(data.flatten(order='F'), config, rms, shape)
	cb = np.reshape(mv, data.shape, order='F')
	return cb


def fellwalker(data not None, config not None, rms):
	data = data.copy()
	shape = np.asarray(data.shape,dtype=np.int32)-1
	mv = _fellwalker(data.flatten(order='F'), config, rms, shape)
	cb = np.reshape(mv, data.shape, order='F')
	return cb
