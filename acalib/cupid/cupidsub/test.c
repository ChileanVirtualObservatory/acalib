#include "cupid.h"
#include "ast.h"
#include "star/hds.h"
#include <stdio.h>
#include <math.h>

CupidPixelSet **cupid_ps_cache = NULL;
int cupid_ps_cache_size = 0;


int main(){
   cupidClumpFind(0,0,NULL,NULL,NULL,NULL,0,NULL,0,0,NULL,NULL,NULL);
   return 0;
}
