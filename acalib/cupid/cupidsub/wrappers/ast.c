#include <stdlib.h>
#include "ast.h"

void *astMalloc_( size_t size, int init, int *status ) {
   // Do not care about security, nor caching... only considering the first argument
   void *result;
   result=malloc(size);
   return result;
}

int *astGetStatusPtr_(){
   // No Status of AST (no AST at all!)
   return NULL;
}

void astAt_( const char *routine, const char *file, int line, int forn,
             int *status) {
   //Do not care for where the problem is (maybe can be used later for logging)
}


void *astFree_( void *ptr, int *status ){
    free(ptr);
    return ptr;
}
