#include "star/hds.h"

/*=================================*/
/* datAlter - Alter size of object */
/*=================================*/

int
datAlter(const HDSLoc    *locator,
         int       ndim,
         const hdsdim    dims[],
         int       *status){
    return -1;
}

/*==========================*/
/* datAnnul - Annul locator */
/*==========================*/

int
datAnnul(HDSLoc    **locator,
         int       *status){
    return -1;
}


/*===========================================*/
/* datCell - Locate a "cell" (array element) */
/*===========================================*/

int
datCell( const HDSLoc   *locator1,
         int      ndim,
         const hdsdim   subs[],
	 HDSLoc   **locator2,
         int      *status){
    return -1;
}


/*================================*/
/* datFind - Find named component */
/*================================*/

int
datFind( const HDSLoc   *locator1,
         const char     *name_str,
	 HDSLoc   **locator2,
         int      *status ){
    return -1;
}

/*============================================*/
/* datNew - Create new component */
/*============================================*/

int
datNew( const HDSLoc    *locator,
        const char      *name_str,
        const char      *type_str,
        int       ndim,
        const hdsdim    dims[],
        int       *status){
    return -1;
}


/*====================================*/
/* datPutD - Write _DOUBLE primitives */
/*====================================*/

int
datPutD( const HDSLoc *locator,
         int       ndim,
         const hdsdim dims[],
         const double    values[],
         int       *status){
    return -1;
}

/*===============================*/
/* datSize - Enquire object size */
/*===============================*/

int
datSize(const HDSLoc *locator,
        size_t *size,
        int *status ){
    return -1;
}


/*===================================*/
/* datTemp - Create temporary object */
/*===================================*/

int
datTemp(const char      *type_str,
        int       ndim,
        const hdsdim    dims[],
        HDSLoc    **locator,
        int       *status){
    return -1;
}

