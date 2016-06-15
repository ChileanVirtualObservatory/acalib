
#include <stdio.h>
#include "mers.h"
#include "sae_par.h"


void msgBlankif( msglev_t prior, int *status ){
  /*  Check the inherited global status. */
  if (*status != SAI__OK) return;

  /* Deliver the message with specified priority */
  msgOutif( prior, "MSG_BLANK", " ", status );

}

void msgOutif( msglev_t prior, const char * param, const char * text, int * status) {
  if (*status != SAI__OK) return;
  // param ignored for the moment
  printf("Priority %d: %s\n",prior,text);
}

void msgSetd( const char *token,
              double dvalue ){
   // Not Implementing an error system TODO: maybe a exception?
}

void msgSetc( const char *token,
              const char *cvalue ){
   // Not Implementing an error system TODO: maybe a exception?
}

void msgSeti( const char *token,
              int ivalue ){
   // Not Implementing an error system TODO: maybe a exception?
}

void errRep( const char *param,
             const char *text,
             int *status ){
    if (*status != SAI__OK) return;
    // param ignored for the moment
    printf("ERROR: %s\n",text);
}

