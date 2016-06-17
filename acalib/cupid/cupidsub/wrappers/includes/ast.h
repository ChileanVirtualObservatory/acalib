#if !defined(AST_INCLUDED)
#define AST_INCLUDED

#include <stdlib.h>
#include <stdarg.h>
#include <float.h>
#include <stdio.h>
#include "Python.h"

#define AST__BAD (-(DBL_MAX))

#define astOK 1

typedef PyObject AstObject;

typedef struct AstFrameSet_t {

} AstFrameSet;

typedef PyDictObject AstKeyMap;


typedef struct AstRegion_t {

} AstRegion;

typedef struct AstFrame_t {

} AstFrame;

typedef struct AstMapping_t {

} AstMapping;

void *astMalloc( size_t size);
void *astFree( void *ptr);
int astMapGet0D( AstKeyMap *map, const char *key, double *value);
int astMapGet0A( AstKeyMap *map, const char *key, AstObject **obj);
void astMapPut0D( AstKeyMap *map, const char *key, double value, const char *comment);
int astMapGet0C( AstKeyMap *map, const char *key, const char **value);
AstObject *astAnnul( AstObject *obj);
void astMapRemove( AstKeyMap *this, const char *key);
void astMapPut0I( AstKeyMap *map, const char *key, int value, const char *comment);
int astMapGet0I( AstKeyMap *map, const char *key, int *value);
void *astGrow( void *ptr, int n, size_t size);
int astSscanf( const char *str, const char *fmt, ...);                            
int astMapSize( AstKeyMap *map);
const char *astMapKey( AstKeyMap *map, int key);

#endif
