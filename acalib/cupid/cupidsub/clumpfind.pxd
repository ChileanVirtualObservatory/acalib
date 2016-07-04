cimport numpy as cnp
#from cpython.cobject cimport PyDictObject
#from cpython cimport PyCDictObject as PyDictObject
#from cpython cimport *
#$from cpython.cobject cimport *
#from cpython.dict cimport PyDictObject

cdef extern from "Python.h":
	ctypedef struct PyObject
	ctypedef struct PyDictObject

cdef extern from "./wrappers/includes/ast.h":
	ctypedef PyDictObject AstKeyMap

cdef extern from "./cupid.h":
	cdef int *cupidClumpFind( int type, int ndim, int *slbnd, int *subnd, 
	void *ipd, double *ipv, double rms, AstKeyMap *config, int velax, 
	int perspectrum, double beamcorr[ 3 ], int *backoff, int *status ) 
