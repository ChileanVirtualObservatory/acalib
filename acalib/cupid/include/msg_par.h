#ifndef MSG_PAR_DEFINED
#define MSG_PAR_DEFINED

typedef enum msglev_t {
  MSG__NONE  = 0, /*   No messages at all */
  MSG__QUIET = 1, /*   Quiet conditional message output level */
  MSG__NORM  = 2, /*   Normal conditional message output level */
  MSG__VERB  = 3, /*   Verbose conditional message output level */
  MSG__DEBUG = 4, /*   Debug conditional message output level */
  MSG__DEBUG1 = 5,/*   Numbered debug levels for flexibility */
  MSG__DEBUG2 = 6,
  MSG__DEBUG3 = 7,
  MSG__DEBUG4 = 8,
  MSG__DEBUG5 = 9,
  MSG__DEBUG6 = 10,
  MSG__DEBUG7 = 11,
  MSG__DEBUG8 = 12,
  MSG__DEBUG9 = 13,
  MSG__DEBUG10 = 14,
  MSG__DEBUG11 = 15,
  MSG__DEBUG12 = 16,
  MSG__DEBUG13 = 17,
  MSG__DEBUG14 = 18,
  MSG__DEBUG15 = 19,
  MSG__DEBUG16 = 20,
  MSG__DEBUG17 = 21,
  MSG__DEBUG18 = 22,
  MSG__DEBUG19 = 23,
  MSG__DEBUG20 = 24,
  MSG__ALL     = 25 /* All messages */
} msglev_t;

#endif  /* MSG_PAR_DEFINED */

