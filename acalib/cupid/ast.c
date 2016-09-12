#include "ast.h"

void *astMalloc( size_t size) {
   // Do not care about security, nor caching... only considering the first argument
   void *result;
   result=malloc(size);
   return result;
}


void *astCalloc( size_t nmemb, size_t size) {
   // Do not care about security, nor caching... only considering the first argument
   void *result;
   result=calloc(nmemb, size);
   return result;
}


void *astFree( void *ptr){
    free(ptr);
    return ptr;
}

AstObject *astAnnul( AstObject *obj){
    //Py_DECREF((PyObject *)obj);
    return NULL;
}

void *astGrow( void *ptr, int n, size_t size){
    ptr=realloc(ptr, size*n);
    return ptr;
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
   PyObject *pkey=PyBytes_FromString(key);
   PyObject *pmap=(PyObject *)map;
   pvalue=PyDict_GetItem(pmap,pkey);
   Py_DECREF(pkey);
   return pvalue;
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
   return 0;
   //*obj=dictGet(map,key);
   //if (obj==NULL) return 0;
   //return 1;
}

void dictPut(AstKeyMap *map, const char *key,PyObject *pvalue){
   PyObject *pkey=PyBytes_FromString(key);
   PyObject *pmap=(PyObject *)map;
   PyDict_SetItem(pmap,pkey,pvalue);
   Py_DECREF(pkey);
}


void astMapPut0D( AstKeyMap *map, const char *key, double value, const char *comment){
   // Comment is neglected 
   PyObject *pvalue=PyFloat_FromDouble(value);
   dictPut(map,key,pvalue);
   Py_DECREF(pvalue);
}

int astMapGet0C( AstKeyMap *map, const char *key, const char **value){
   PyObject *pvalue;
   pvalue=dictGet(map,key);
   if (pvalue==NULL) return 0;
   *value=PyBytes_AsString(pvalue);
   Py_DECREF(pvalue);
   return 1;
}


void astMapRemove(AstKeyMap *map, const char *key){
   PyObject *pkey=PyBytes_FromString(key);
   PyObject *pmap=(PyObject *)map;
   PyDict_DelItem(pmap,pkey);
   Py_DECREF(pkey);
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
    PyObject *pvalue=PyLong_FromLong(value);
    dictPut(map,key,pvalue);
    Py_DECREF(pvalue);
}


int astMapSize( AstKeyMap *map){
    return map->ma_used;
}

const char *astMapKey( AstKeyMap *map, int key){
    PyObject *pmap=(PyObject *)map;
    PyObject *plist=PyDict_Keys(pmap);
    PyObject *pitem=PyList_GetItem(plist,key);
    Py_DECREF(plist);
    const char *retval=PyBytes_AsString(pitem);
    Py_DECREF(pitem);
    return retval;
}


