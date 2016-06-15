#ifndef KAPLIBS_DEFINED
#define KAPLIBS_DEFINED
#include <stdio.h>
#include "star/grp.h"
#include "ast.h"


void kpg1Rgndf( const char *, size_t, size_t, const char *, Grp **, size_t *, int * );
void kpg1Gtgrp( const char *, Grp **, size_t*, int * );
void kpg1Kygrp( AstKeyMap *, Grp **, int * );
void kpg1Asget( int, int, int, int, int, int *, int *, int *, AstFrameSet **, int * );
void kpg1Kymap( const Grp *, AstKeyMap **, int * );


#endif  /* CVG_DEFINED */

