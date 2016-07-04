/*
*+
*  Name:
*     cupid.h

*  Purpose:
*     Standard include file for CUPID

*  Language:
*     {routine_language}

*  Copyright:
*     Copyright (C) 2005-2006 Particle Physics & Astronomy Research Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software; you can redistribute it and/or
*     modify it under the terms of the GNU General Public License as
*     published by the Free Software Foundation; either version 2 of
*     the License, or (at your option) any later version.
*
*     This program is distributed in the hope that it will be
*     useful, but WITHOUT ANY WARRANTY; without even the implied
*     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
*     PURPOSE. See the GNU General Public License for more details.
*
*     You should have received a copy of the GNU General Public License
*     along with this program; if not, write to the Free Software
*     Foundation, Inc., 51 Franklin Street,Fifth Floor, Boston, MA
*     02110-1301, USA

*  Authors:
*     DSB: David S. Berry (UCLan)
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  History:
*     14-JAN-2009 (TIMJ):
*        Remove ilevel from APIs
*     {enter_changes_here}

*  Bugs:
*     {note_any_bugs_here}

*-
*/

#if !defined( CUPID_INCLUDED )   /* Include this file only once */
#define CUPID_INCLUDED

#include "ast.h"
#include "star/grp.h"
#include "star/hds.h"
#include "msg_par.h"

/* Constants */
/* --------- */
#define CUPID__DOUBLE 1
#define CUPID__FLOAT  2

#define CUPID__KBACK  0
#define CUPID__KEDGE  1
#define CUPID__KPEAK  2

#define CUPID__GCNP1  4  /* No. of free parameters in a 1D GaussClump clump */
#define CUPID__GCNP2  7  /* No. of free parameters in a 2D GaussClump clump */
#define CUPID__GCNP3 11  /* No. of free parameters in a 3D GaussClump clump */

#define CUPID__CFNULL -1 /* Unassigned pixel flag */

#define CUPID__CONFIG  "NOALG_CONFIG" /* Key for config params which have no
                                         algorithm name */

/* Type Definitions */
/* ---------------- */

/* A structure used to describe an iterator that iterates round the
   pixels in a rectangular subsection of a 3D array. */
typedef struct CupidBoxIter {
   int done;          /* Have all required pixels been visited? */
   int xx0;           /* X grid index of current array element */
   int xx1;           /* Y grid index of current array element */
   int xx2;           /* Z grid index of current array element */
   int i;             /* Index of current array element in overlap */
   int i2;            /* Index of current xy plane in overlap */
   int i3;            /* Index of current x row in overlap */
   int lbnd0;         /* Lower X bounds of rectangle being processed */
   int lbnd1;         /* Lower Y bounds of rectangle being processed */
   int lbnd2;         /* Lower Z bounds of rectangle being processed */
   int ubnd0;         /* Upper X bounds of rectangle being processed */
   int ubnd1;         /* Upper Y bounds of rectangle being processed */
   int ubnd2;         /* Upper Z bounds of rectangle being processed */
   int xsize;         /* The X dimension of the array */
   int xysize;        /* No. of pixels in an XY plane of the array */
   int gap;           /* Gap between visited pixels */
} CupidBoxIter;


/* A structure holding the global parameters of the GaussClump algorithm
   needed by service functions. */
typedef struct CupidGC  {
   double *data;           /* Pointer to copy of data section being fitted */
   double *res;            /* Pointer to array to receive scale residuals */
   double *resu;           /* Pointer to array to receive unscale residuals */
   double *weight;         /* Pointer to weights for section being fitted */
   double beam_sq;         /* Square of spatial beam FWHM in pixels */
   double lbnd[ 3 ];       /* Lower grid bounds of section being fitted */
   double maxwf;           /* Maximum factor for modifying weights */
   double minwf;           /* Minimum factor for modifying weights */
   double s0p1;            /* Chi-square stiffness parameter s0, minus 1.0 */
   double sa;              /* Chi-square stiffness parameter sa */
   double sb;              /* Chi-square stiffness parameter sb */
   double sc4;             /* Four times chi-square stiffness parameter sc */
   double ubnd[ 3 ];       /* Upper grid bounds of section being fitted */
   double velres_sq;       /* Square of velocity resolution in pixels */
   double wsum;            /* Sum of values in "weight" array */
   double x_max[ 3 ];      /* Grid coords of "ymax" value */
   double ymax;            /* Largest data value in section being fitted */
   int dax[ 3 ];           /* External axis no. indexed by algorithm axis no. */
   int fixback;            /* Is the background to be kept fixed during fit? */
   int maxnf;              /* Max invocation count for calcf */
   int ndim;               /* Number of pixel axes in the data array */
   int nel;                /* Number of pixels in section being fitted */
   int nf;                 /* The invocation count from calcf */
   int npar;               /* No of free params in fit (inc. background level) */
   int nwf;                /* Number of times to modifiy the weight array */

   double initpars[ CUPID__GCNP3 ];
   double pars[ CUPID__GCNP3 ];
   double *initmodel;
   double *model;
   double *resids;
   double chisq;
   int slbnd[3];           /* Lower pixel bounds of user supplied NDF */

} CupidGC;

/* Structure used to describe a set of contiguous pixels in the ClumpFind
   and FellWalker algorithms. */
typedef struct CupidPixelSet {
   int index;          /* The index used to identify pixels in the set */
   int lbnd[ 3 ];      /* Lower GRID bounds of the set bounding box */
   int ubnd[ 3 ];      /* Upper GRID bounds of the set bounding box */
   double vpeak;       /* Peak pixel value */
   int peak[ 3 ];      /* Peak pixel GRID coords */
   int pop;            /* Number of pixels assigned to the pixel set */
   int edge;           /* Does the PixelSet touch an edge of the data array? */
   int *nebs;          /* List of indices of neighbouring clumps */
   int nneb;           /* Length of "nebs" list */
   double *cols;       /* FW: data value at the "col" between two neigbours */
   int lneb;           /* Cached clump index */
   int lnebi;          /* Index of cached clump index within "nebs" and "cols" */
} CupidPixelSet;


/* A structure used to store information about a group of clumps (used by
   CLUMPINFO). */
typedef struct CupidClumpInfo {
   int init;           /* Has the structure been initialised? */
   AstFrameSet *iwcs;  /* WCS FrameSet from main NDF */
   int lbnd[ 3 ];      /* Lower pixel index bounds of bounding box */
   int ubnd[ 3 ];      /* Upper pixel index bounds of bounding box */
   int npix;           /* Number of pixel axes */
   int nwcs;           /* Number of WCS axes */
} CupidClumpInfo;

/* A structure used to store information required by cupidFindback0 (the
   thread worker function used by the FINDBACK command). */
typedef struct CupidFindback0Data {
   int islice;         /* Slice index */
   int nslice;         /* Number of slices to process */
   int type;           /* Integer identifier for data type */
   int ndim;           /* Total number of pixel axes in NDF */
   int box[ 3 ];       /* Dimensions of each cell in pixels */
   double rms;         /* Global rms error in data */
   void *ipd1;         /* Pointer to input Data array */
   void *ipd2;         /* Pointer to output Data array */
   int slice_dim[ 3 ]; /* Dimensions of each significant slice axis */
   int slice_lbnd[ 3 ];/* Lower bounds of each significant slice axis */
   int newalg;         /* Use experimental algorithm variations? */
   int slice_size;     /* Number of pixels in each slice */
   float wlim;         /* Min. frac. of good pixels in a filter box */
} CupidFindback0Data;



/* Function macros */
/* --------------- */

/* A macro to return the size of a cupid data type. */
#define cupidSize( type, fun ) \
   ( ( type == CUPID__DOUBLE ) ? \
      sizeof( double ) \
\
   : ( ( type == CUPID__FLOAT ) ? \
       sizeof( float ) \
\
   : ( ( *status == SAI__OK ) ? \
       *status = SAI__ERROR, \
       msgSeti( "TYPE", type ), \
       errRep( "CUPIDSIZE_ERR1", fun ": Invalid \"type\" " \
               "value (^TYPE) supplied (CUPID programming error).", \
               status ), 0 : 0 ) ) )




#define cupidTestBnd \
{int itea, iteax, iteay, iteaz; \
   itea = -1; \
   for( iteaz = 1; iteaz <= dims[2]; iteaz++ ) { \
      for( iteay = 1; iteay <= dims[1]; iteay++ ) { \
         for( iteax = 1; iteax <= dims[0]; iteax++ ) { \
            itea++; \
            if( cupidMergeSet( ipa[ itea ] ) ) == ps->index ) { \
               if( iteax < ps->lbnd[ 0 ] || iteax > ps->ubnd[ 0 ] ||  \
                   iteay < ps->lbnd[ 1 ] || iteay > ps->ubnd[ 1 ] ||  \
                   iteaz < ps->lbnd[ 2 ] || iteaz > ps->ubnd[ 2 ] ) { \
                  printf("Pixel %d [%d %d %d] has ipa %d but is not in " \
                         "bounding box (%d:%d,%d:%d)\n", itea, iteax, iteay,  \
                         iteaz, ipa[itea], ps->lbnd[ 0 ], ps->ubnd[ 0 ], \
                         ps->lbnd[ 1 ], ps->ubnd[ 1 ]); \
               } \
            } \
         } \
      } \
   }





/* A set of macros which invoke the corresponding AST memory management
   functions, but which specify the size of a data element using a cupid
   data type rather than a C data type. The "fun" macro parameter should
   be a quoted string holding the name of the function to be included in
   any error messages. */

#define cupidStore( mem, ptr, nel, type, fun ) \
   astStore( mem, ptr, nel*cupidSize( type, fun ) );


/* PixelSet cache used by the ClumpFind algorithm. */
/* ----------------------------------------------- */
/* Pointer to an array holding a list of PixelSet pointers. The PixelSet
   structures in this list have been created previously by cupidCFMakePS
   but have subsequently been freed (using cupidCFFreePS) and are
   currently not being used for anything, and so can be re-issued by
   cupidCFMakePS, thus avoiding the overhead of frequenct memory
   allocation. */
extern CupidPixelSet **cupid_ps_cache;

/* This is the length of the cupid_ps_cache array. */
extern int cupid_ps_cache_size;

/* Function prototypes */
/* ------------------  */
AstKeyMap *cupidRetrieveConfig( HDSLoc *, int * );
AstRegion *cupidEllipseDesc( AstFrame *, int[ 2 ], double[ 3 ], double, double, double, double, double, double, float[ 4 ], int, int *, AstMapping *, AstFrame *, AstMapping *, int * );
AstRegion *cupidPolygonDesc( double *ipd, int velax, double *peak, int space_axes[ 2 ], int ndim, int *lbnd, int *ubnd, AstMapping *wcsmap, AstFrame *space_frm, AstMapping *space_map, int *status );
CupidBoxIter *cupidBoxIterator( CupidBoxIter *, int[3], int[3], int[3], int, int * );
CupidPixelSet *cupidCFDeletePS( CupidPixelSet *, int * );
CupidPixelSet *cupidCFFreePS( CupidPixelSet *, int *, int, int * );
CupidPixelSet *cupidCFMakePS( int, int * );
HDSLoc *cupidClumpFind( int, int, int *, int *, void *, double *, double, AstKeyMap *, int, int, double[3], int *, int * );
HDSLoc *cupidFellWalker( int, int, int *, int *, void *, double *, double, AstKeyMap *, int, int, double[3], int * );
HDSLoc *cupidGaussClumps( int, int, int *, int *, void *, double *, double, AstKeyMap *, int, double[3], int * );
HDSLoc *cupidReinhold( int, int, int *, int *, void *, double *, double, AstKeyMap *, int, double[3], int * );
double *cupidCFLevels( AstKeyMap *, double, double, double, int *, int * );
double *cupidClumpDesc( int, int, AstMapping *, AstFrame *, const char *, double[ 3 ], int, int, int, double *, const char ***, const char ***, int *, int *, char **, AstRegion **, int * );
double cupidConfigD( AstKeyMap *, const char *, double, int * );
double cupidConfigRMS( AstKeyMap *, const char *, double, double, int * );
double cupidGCChiSq( int, double *, int, int, int * );
double cupidGCModel( int, double *, double *, int, int, int, int * );
float cupidRanVal( int, float[2], int * );
int *cupidRCA( int *, int *, int, int[ 3 ], int[ 3 ], double, int, int, int, int, int * );
int *cupidRCA2( int *, int *, int, int[ 3 ], int[ 3 ], int * );
int cupidCFErode( CupidPixelSet *, int *, int, int *, int[3], int, int, CupidPixelSet **, int * );
int cupidCFXtend( CupidPixelSet *, CupidPixelSet *, int *, int, int *, int[3], int, CupidPixelSet **, int * );
int cupidConfigI( AstKeyMap *, const char *, int, int * );
int cupidDefMinPix( int, double *, double, double, int * );
int cupidNextIt( CupidBoxIter *, int[3], int *, int * );
int cupidRFillClumps( int *, int *, int, int, int[ 3 ], int[ 3 ], int, int * );
void cupidCFAddPixel( int *, CupidPixelSet *, int, int[3], double, int, int * );
void cupidCFIdl( CupidPixelSet *, int *, int, int *, int[3], int, CupidPixelSet **, int * );
void cupidCFMerge( CupidPixelSet *, CupidPixelSet *, int *, int[3], int *, int **, int, CupidPixelSet **, int * );
void cupidCFNebs( int *, int, int x[], int, int [3], int[3], int, int, int, int *, int *, int[27], int *, int *, CupidPixelSet **, int * );
void cupidCFXfer( CupidPixelSet *, CupidPixelSet *, int *, int[3], int * );
void cupidClumpInfo1( HDSLoc *, CupidClumpInfo *, int * );
void cupidDatCopy( HDSLoc *, HDSLoc *, int * );
void cupidDumpD( double *, int, int *, int *, const char *, int * );
void cupidDumpF( float *, int, int *, int *, const char *, int * );
void cupidDumpI( int *, int, int *, int *, const char *, int * );
void cupidEdges( float *, int, int[3], int[3], float, float, int * );
void cupidFindback0( void *, int * );
void cupidGCListClump( int, int, double *, double, int *, double, int * );
void cupidGCNdfClump( HDSLoc **, double, double *, double, int, int *, int *, int, double *, int *, int *, int, int *, AstKeyMap *, int, int * );
void cupidGCcalcf( int, double *, int *, double * );
void cupidGCcalcg( int, double *, int *, double * );
void cupidREdges( int, double *, int *, int *, int, double, double, double, double, int * );
void cupidRFillLine( int *, int *, int, int, int[ 3 ], int[ 3 ], int[ 3 ], int, int, int, int, int *[3], int * );
void cupidStoreClumps( const char *, const char *, int, HDSLoc *, HDSLoc *, int, int, int, int, int, double[ 3 ], const char *, int, AstFrameSet *, const char *, Grp *, FILE *, int *, int * );
void cupidStoreConfig( HDSLoc *, AstKeyMap *, int * );

void findback( int * );
void findclumps( int * );
void cupidhelp( int * );
void makeclumps( int * );
void extractclumps( int * );
void clumpinfo( int * );

void cupid_mon( int * );
