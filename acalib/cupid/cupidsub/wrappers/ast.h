#if !defined(AST_INCLUDED)
#define AST_INCLUDED

#include <stdlib.h>

typedef struct AstObject_t{

} AstObject;

typedef struct AstFrameSet_t {

} AstFrameSet;

typedef struct AstKeyMap_t {

} AstKeyMap;


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
int astMapGet0I( AstKeyMap *map, const char *key, int *value);



#endif
