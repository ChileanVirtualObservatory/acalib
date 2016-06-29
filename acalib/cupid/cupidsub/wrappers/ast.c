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

PyObject *dictGet(AstKeyMap *map, const char *key){
   PyObject *pvalue;
   PyString *pkey=PyString_FromString(key);
   pvalue=PyDict_GetItem(map,pkey);
   Py_DECREF(pkey);
   return pvalue
}

int astMapGet0D( AstKeyMap *map, const char *key, double *value){
   PyObject *pvalue;
   pvalue=dictGet(map,key);
   if (pvalue==NULL) return 0;
   *value=PyFloat_AsDouble(pvalue);
   Py_DECREF(pvalue);
   return 1;
}

int astMapGet0A( AstKeyMap *map, const char *key, AstObject **obj){
   *obj=dictGet(map,key);
   if (obj==NULL) return 0;
   return 1;
}

void astMapPut0D( AstKeyMap *map, const char *key, double value, const char *comment){
   
}

int astMapGet0C( AstKeyMap *map, const char *key, const char **value){
   PyObject *pvalue;
   pvalue=dictGet(map,key);
   if (pvalue==NULL) return 0;
   *value=PyUnicode_AsUTF8(pvalue);
   Py_DECREF(pvalue);
   return 1;
}


void astMapRemove( AstKeyMap *this, const char *key){

}

int astMapGet0I(AstKeyMap *map, const char *key, int *value){
   PyObject *pvalue;
   pvalue=dictGet(map,key);
   if (pvalue==NULL) return 0;
   *value=PyLong_AsLong(pvalue);
   Py_DECREF(pvalue);
   return 1;
}

void astMapPut0I( AstKeyMap *map, const char *key, int value, const char *comment){
   return 1;
}


int astMapSize( AstKeyMap *map){
    return -1;

}

const char *astMapKey( AstKeyMap *map, int key){
    return NULL;
}


