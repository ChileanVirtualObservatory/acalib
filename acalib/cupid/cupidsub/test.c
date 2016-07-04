#include "cupid.h"
#include "ast.h"
#include "star/hds.h"
#include "Python.h"
#include <stdio.h>
#include <math.h>

CupidPixelSet **cupid_ps_cache = NULL;
int cupid_ps_cache_size = 0;

int main(){
  
   Py_Initialize();
   cupidClumpFind(0,0,NULL,NULL,NULL,NULL,0,NULL,0,0,NULL,NULL,NULL);
   PyObject *pmap=PyDict_New();
   AstKeyMap *map=(AstKeyMap *)pmap;
   astMapPut0D(map,"ONE", 1.0, NULL);
   printf("%d\n",astMapSize(map));
   int ival;
   double dval;
   astMapGet0D(map,"ONE",&dval);
   printf("%f\n",dval);
   astMapPut0I(map,"TWO", 2, NULL);
   printf("%d\n",astMapSize(map));
   astMapGet0I(map,"TWO",&ival);
   printf("%i\n",ival);
   return 0;
}
