#ifndef MERS_DEFINED
#define MERS_DEFINED

#include "msg_par.h"

void msgSetd( const char *token,
              double dvalue );

void msgSeti( const char *token,
              int ivalue );

void msgSetc( const char *token,
              const char *cvalue );

void msgOutif( msglev_t prior,
               const char *param,
               const char *text,
               int *status );

void msgOut( const char *param,
             const char *text,
             int *status );

void msgOutf( const char *param,
              const char *text,
              int *status,
              ... ) __attribute__((format (printf, 2, 4 )));

void msgBlankif( msglev_t prior, int *status );

void msgSetr( const char *token, float rvalue );

int msgFlevok( msglev_t filter, int *status );

void errRep( const char *param,
             const char *text,
             int *status );

void errAnnul( int *status );

void errRepf( const char *param,
              const char *text,
              int *status,
              ... ) __attribute__((format (printf, 2, 4 )));

void errBegin( int *status );

void errEnd( int *status );



#endif  /* MERS_DEFINED */

