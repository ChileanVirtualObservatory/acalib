cimport numpy as cnp
#from cpython.cobject cimport PyDictObject
#from cpython cimport PyCDictObject as PyDictObject
#from cpython cimport *
#$from cpython.cobject cimport *
#from cpython.dict cimport PyDictObject

cdef extern from "Python.h":
	cdef struct PyObject
	cdef struct PyDictObject

cdef extern from "./cupidsub/wrappers/includes/ast.h":
	cdef struct AstKeyMap

cdef extern from "./cupidsub/cupid.h":
	cdef int *cupidClumpFind( int type, int ndim, int *slbnd, int *subnd, 
	void *ipd, double *ipv, double rms, AstKeyMap *config, int velax, 
	int perspectrum, double beamcorr[ 3 ], int *backoff, int *status )
