#ifndef GRP_DEFINED
#define GRP_DEFINED

/* Note that GRP__NOIDs role in C is to use a NULL pointer. */

/* Maximum length of a group expression. */
enum { GRP__SZGEX  = 255 };

/* Length of a name within a group. */
enum { GRP__SZNAM  = 255 };

/* Max. length of a group type */
enum { GRP__SZTYP  = 80 };

/* Max. length of a file name. */
enum { GRP__SZFNM  = 256 };

/* Max. number of groups which can be used simultaneously. */
enum { GRP__MAXG  = 2048 };


typedef struct Grp_t{

} Grp;

#endif  /* GRP_DEFINED */

