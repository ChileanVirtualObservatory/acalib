#ifndef CVG_DEFINED
#define CVG_DEFINED

#include "fitsio.h"

void cvgAssoc( const char *param, const char *mode, fitsfile **fptr, int *blockf, int *status );
void cvgWhisr( int ndf, fitsfile *fptr, int *status );
void cvgClose( fitsfile **fptr, int *status );


#endif  /* CVG_DEFINED */

