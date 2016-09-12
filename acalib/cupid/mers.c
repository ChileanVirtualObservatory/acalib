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

void msgOut( const char * param, const char * text, int * status ) {
  if (*status != SAI__OK) return;  
}

void msgOutf( const char * param, const char * text, int *status, ... ) {
  if (*status != SAI__OK) return;
}

void msgSetd( const char *token,
              double dvalue ){
   //TODO: implement token-based print
}

void msgSetc( const char *token,
              const char *cvalue ){
   //TODO: implement token-based print
}

void msgSeti( const char *token,
              int ivalue ){
   //TODO: implement token-based print
}

void msgSetr( const char *token, float rvalue ){
  //TODO: implement token-based print
}

int msgFlevok( msglev_t filter, int *status ){
  //TODO: implement token-based print
}

void errRep( const char *param,
             const char *text,
             int *status ){
    if (*status != SAI__OK) return;
    // param ignored for the moment
    printf("ERROR: %s\n",text);
}

void errAnnul( int * status ){

}

