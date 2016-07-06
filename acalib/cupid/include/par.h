#ifndef PAR_DEFINED
#define PAR_DEFINED

/* Null parameter value */
enum { PAR__NULL            	= 146703163 };	/* messid=103 */

void parGet0c( const char *param,
               char *value,
               int value_length,
               int *status );
void parDef0l( const char *param,
               int value,
               int *status );
void parGet0l( const char *param,
               int *value,
               int *status );
void parChoic( const char *param,
               const char *defaul,
               const char *opts,
               int null,
               char *value,
               int value_length,
               int *status );
void parDef0d( const char *param,
               double value,
               int *status );
void parGet0d( const char *param,
               double *value,
               int *status );
void parPut0i( const char *param,
               int value,
               int *status );

#endif  /* PAR_DEFINED */

