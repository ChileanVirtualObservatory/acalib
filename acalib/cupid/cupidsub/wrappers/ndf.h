#ifndef NDF_DEFINED
#define NDF_DEFINED

/*  Maximum number of NDF dimensions.                                       */
#define NDF__MXDIM 7
#include "star/hds.h"

void ndfBegin( void );
void ndfCget( int indf,
              const char *comp,
              char *value,
              int value_length,
              int *status );
void ndfDim( int indf,
             int ndimx,
             int dim[],
             int *ndim,
             int *status );
void ndfMsg( const char *token,
             int indf );
void ndfMtype( const char *typlst,
               int indf1,
               int indf2,
               const char *comp,
               char *itype,
               int itype_length,
               char *dtype,
               int dtype_length,
               int *status );
void ndfMap( int indf,
             const char *comp,
             const char *type,
             const char *mmod,
             void *pntr[],
             int *el,
             int *status );
void ndfState( int indf,
               const char *comp,
               int *state,
               int *status );
void ndfProp( int indf1,
              const char *clist,
              const char *param,
              int *indf2,
              int *status );
void ndfSbad( int bad,
              int indf,
              const char *comp,
              int *status );
void ndfStype( const char *ftype,
               int indf,
               const char *comp,
               int *status );
void ndfXnew( int indf,
              const char *xname,
              const char *type,
              int ndim,
              const int dim[],
              HDSLoc **loc,
              int *status );
void ndfHdef( int indf,
              const char *appn,
              int *status );
void ndfLoc( int indf,
             const char *mode,
             HDSLoc ** loc,
             int *status );
void ndfXstat( int indf,
               const char *xname,
               int *there,
               int *status );
void ndfXdel( int indf,
              const char *xname,
              int *status );
void ndfAnnul( int *indf,
               int *status );
void ndfEnd( int *status );
void ndfDelet( int *indf,
               int *status );

#endif  /* NDF_DEFINED */

