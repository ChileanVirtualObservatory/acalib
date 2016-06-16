#include "ast.h"

void *astMalloc( size_t size) {
   // Do not care about security, nor caching... only considering the first argument
   void *result;
   result=malloc(size);
   return result;
}


void *astFree( void *ptr){
    free(ptr);
    return ptr;
}

AstObject *astAnnul( AstObject *obj){
   return NULL;
}

void *astGrow( void *ptr, int n, size_t size){
    return NULL;
}


int astSscanf( const char *str, const char *fmt, ...){
   /* Initialise the variable argument list pointer. */
   va_list args;  
   va_start( args, fmt );
   return vsscanf(str,fmt,args);
}

// KEYMAP

int astMapGet0D( AstKeyMap *map, const char *key, double *value){
   // Not implemented
   return -1;
}

int astMapGet0A( AstKeyMap *map, const char *key, AstObject **obj){
   // Not implemented
   return -1;
}

void astMapPut0D( AstKeyMap *map, const char *key, double value, const char *comment){
   
}

int astMapGet0C( AstKeyMap *map, const char *key, const char **value){
    return -1;
}


void astMapRemove( AstKeyMap *this, const char *key){

}

int astMapGet0I(AstKeyMap *map, const char *key, int *value){
    return -1;
}

void astMapPut0I( AstKeyMap *map, const char *key, int value, const char *comment){
    
}


int astMapSize( AstKeyMap *map){
    return -1;

}

const char *astMapKey( AstKeyMap *map, int key){
    return NULL;
}


