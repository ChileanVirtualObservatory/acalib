import numpy as np
cimport cupid
from libc.stdio cimport printf

cdef int[:] gc(double[::1] data, double[::1] variance,
			   int[::1] shape, cupid.AstKeyMap *config, 
			   double rms, int velax):
	cdef:
		int[::1] _slbnd = np.zeros(shape.size, dtype=np.int32)

	cdef:
		int dtype = cupid.CUPID__DOUBLE
		int ndim = shape.size
		int *slbnd = &_slbnd[0]
		int *subnd = <int *> &shape[0]
		void *ipd = &data[0]
		double *ipv = &variance[0]

	cdef:
		cupid.HDSLoc *ndfs
		double beamcorr[3]
		int status = 0
		int[:] output = np.zeros(shape.size, dtype=np.int32) #dummy

	#ndfs = cupid.cupidGaussClumps(dtype, ndim, slbnd, subnd,
	#							  ipd, ipv, rms, config, velax,
	#							  beamcorr, &status)

	#somefunction to get output array.

	return output

def gaussclumps(data, variance, config, rms, velax):
	cdef:
		cupid.AstKeyMap *aconfig = cupid.astKeyMap(" ")

	cdef:
		# Variables used to test the ast library.
		int int_value
		double double_value
		const char *string_value

	for key, value in config.items():
		key = bytes(key, "ascii")
		if type(value) is int:
			cupid.astMapPut0I(aconfig, key, value, NULL)
		elif type(value) is float:
			cupid.astMapPut0D(aconfig, key, value, NULL)
		elif type(value) is str:
			cupid.astMapPut0C(aconfig, key, bytes(value, "ascii"), NULL)
		else:
			raise ValueError("Value for key " + repr(key.decode()) 
							 + " should be of type int, float or str")

	shape = np.asarray(data.shape, dtype=np.int32) - 1
	data = data.flatten(order='F')
	variance = variance.flatten(order='F')

	# Test AstKeyMap
	for key, value in config.items():
		key = bytes(key, "ascii")
		printf("%s: ", <char *> key)
		if type(value) is int:
			cupid.astMapGet0I(aconfig, key, &int_value)
			printf("%d\n", int_value)
		elif type(value) is float:
			cupid.astMapGet0D(aconfig, key, &double_value)
			printf("%g\n", double_value)
		elif type(value) is str:
			cupid.astMapGet0C(aconfig, key, &string_value)
			printf("%s\n", string_value)

	clumps = gc(data, variance, shape, aconfig, rms, velax)