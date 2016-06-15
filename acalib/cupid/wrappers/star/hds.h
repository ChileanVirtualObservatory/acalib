#ifndef HDS_DEFINED
#define HDS_DEFINED

typedef struct HDSLoc_t {

} HDSLoc;

int
datFind( const HDSLoc   *locator1,
         const char     *name_str,
    HDSLoc   **locator2,
         int      *status );
/*=======================*/
/* datCopy - copy object */
/*=======================*/

int
datCopy(const HDSLoc    *locator1,
   const HDSLoc    *locator2,
        const char      *name_c,
        int       *status );
/*==========================*/
/* datAnnul - Annul locator */
/*==========================*/

int
datAnnul(HDSLoc    **locator,
         int       *status);


#endif  /* HDS_DEFINED */

