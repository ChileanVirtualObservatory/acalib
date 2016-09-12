cimport numpy as cnp
#from cpython.cobject cimport PyDictObject
#from cpython cimport PyCDictObject as PyDictObject
#from cpython cimport *
#$from cpython.cobject cimport *
#from cpython.dict cimport PyDictObject

ctypedef PyDictObject AstKeyMap;

cdef extern from "Python.h":
	cdef struct PyObject
	cdef struct PyDictObject

cdef extern from "./cupid.h":
	cdef int *cupidClumpFind( int type, int ndim, int *slbnd, int *subnd,
	void *ipd, double *ipv, double rms, AstKeyMap *config, int velax,
	int perspectrum, double beamcorr[ 3 ], int *backoff, int *status )

	cdef int *cupidFellWalker( int type, int ndim, int *slbnd, int *subnd,
	void *ipd, double *ipv, double rms, AstKeyMap *config, int velax,
    int perspectrum, double beamcorr[ 3 ], int *status )
