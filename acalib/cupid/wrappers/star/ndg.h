#ifndef NDG_DEFINED
#define NDG_DEFINED

#include "star/grp.h"

void ndgNdfas( const Grp *igrp, size_t index, const char mode[], int *indf, int *status );
void ndgHltgh( int new, int *old, int *status );
void ndgHltpv( int new, int *old, int *status );


#endif  /* NDG_DEFINED */

