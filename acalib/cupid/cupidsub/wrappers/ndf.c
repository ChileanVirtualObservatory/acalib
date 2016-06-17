#include "ndf.h"

void ndfBegin( void ){
}
void ndfPlace( const HDSLoc * loc,
               const char *name,
               int *place,
               int *status ){
}
void ndfNew( const char *ftype,
             int ndim,
             const int lbnd[],
             const int ubnd[],
             int *place,
             int *indf,
             int *status ){
}
void ndfMap( int indf,
             const char *comp,
             const char *type,
             const char *mmod,
             void *pntr[],
             int *el,
             int *status ){
}
void ndfUnmap( int indf,
               const char *comp,
               int *status ){
}
void ndfXnew( int indf,
              const char *xname,
              const char *type,
              int ndim,
              const int dim[],
              HDSLoc **loc,
              int *status ){
}
void ndfEnd( int *status ){
}
void ndfCput( const char *value,
              int indf,
              const char *comp,
              int *status ){

}

