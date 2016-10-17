/* mem.h.  Generated from mem.h.in by configure.  */
#if !defined ( STAR_MEM_INCLUDED )   /* Include file only once */
#define STAR_MEM_INCLUDED

/*
*  Name:
*     star/mem.h

*  Purpose:
*     Starlink memory management wrapper routines

*  Description:
*     This include file declares the public interface to the
*     starmem library.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)

*  History:
*     08-FEB-2006 (TIMJ):
*        Original version.
*     03-MAY-2007 (TIMJ):
*        Do not include gc.h unless libgc also available

*  Copyright:
*     Copyright (C) 2006 Particle Physics and Astronomy Research Council.
*     Copyright (C) 2007 Science and Technology Facilities Council.
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
*     You should have received a copy of the GNU General Public
*     License along with this program; if not, write to the Free
*     Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*     MA 02110-1301, USA

*/

/* Make sure we have standard definitions such as size_t */
#include <stdlib.h>

/* Private init routine called by starMemInit macro */
void starMemInitPrivate( int gc_initialized );

/* Initialisation routine. Required to be called from the main
   program (not library). Must be called before any malloc routine
   if the garbage collector is required. This is a #define not a function
   and contains (optionally) the GC_INIT call. We provide a starMemInit()
   function in a .o file (not part of the library) for systems that
   need to link from Fortran.
 */

/* HAVE_GC_H and HAVE_LIBGC are substituted below at configure time */
/* #undef HAVE_LIBGC */
/* #undef HAVE_GC_H */

#if HAVE_GC_H && HAVE_LIBGC
#include <gc.h>
#endif

/* We always call GC_INIT if it is defined, even if we have
   disabled Garbage Collection. This is because this include file
   can not know whether GC is required and it never hurts.
*/

#if HAVE_LIBGC && defined( GC_INIT )
#define MYINIT() GC_INIT()
#define starMemInit() { MYINIT(); starMemInitPrivate(1); }
#else
#define starMemInit() { starMemInitPrivate(0); }
#endif


/* Standard system replacements for malloc/free */


int    starMemIsInitialised( void );
void * starMalloc( size_t size );
void * starMallocAtomic( size_t size );
void * starCalloc( size_t count, size_t size );
void * starRealloc( void * ptr, size_t size );
void   starFree( void * ptr );
void   starFreeForce( void * ptr );

/* STAR_MEM1_INCLUDED */
#endif
