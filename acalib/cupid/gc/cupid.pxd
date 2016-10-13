cdef extern from "ast.h":
	ctypedef struct AstKeyMap:
		pass
	AstKeyMap* astKeyMap(const char *options, ...)
	void astMapPut0D(AstKeyMap *this, const char *key, double value, const char *comment)
	void astMapPut0I(AstKeyMap *this, const char *key, int value, const char *comment)
	void astMapPut0C(AstKeyMap *this, const char *key, const char *value, const char *comment)
	int astMapGet0I(AstKeyMap *this, const char *key, int *value)
	int astMapGet0D(AstKeyMap *this, const char *key, double *value)
	int astMapGet0C(AstKeyMap *this, const char *key, const char **value)
	
cdef extern from "star/hds_types.h":
	ctypedef struct HDSLoc:
		pass

cdef extern from "cupid.h":
	int CUPID__DOUBLE
	HDSLoc* cupidGaussClumps(int type, int ndim, int *slbnd, 
							 int *subnd, void *ipd, double *ipv,
							 double rms, AstKeyMap *config, int velax,
							 double beamcorr[3], int *status)