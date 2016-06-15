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

AstObject *astAnnul( AstObject *obj){
   return NULL;
}

void astMapRemove( AstKeyMap *this, const char *key){

}

int astMapGet0I(AstKeyMap *map, const char *key, int *value){
    return -1;
}

void astMapPut0I( AstKeyMap *map, const char *key, int value, const char *comment){
   
}

//int *astGetStatusPtr_(){
   // No Status of AST (no AST at all!)
//   return NULL;
//}

//void astAt_( const char *routine, const char *file, int line, int forn,
//             int *status) {
   //Do not care for where the problem is (maybe can be used later for logging)
//}
//
//
//AstObject *astMakePointer_( AstObject *this_id, int *status ) {
//   return this_id;
//}

//AstObject *astCheckLock( AstObject *this) {
//   return this;
//}
