#ifndef IRQ_DEFINED
#define IRQ_DEFINED

typedef struct IRQLocs_t {

} IRQLocs;

void irqDelet( int, int * );
void irqNew( int, const char *, IRQLocs **, int * );
void irqAddqn( const IRQLocs *, const char *, int, const char *, int * );
void irqSetqm( const IRQLocs *, int, const char *, int, float *, int *, int * );
void irqRlse( IRQLocs **, int * );


#endif  /* IRQ_DEFINED */

