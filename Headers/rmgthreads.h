#ifndef RMG_THREADS_H
#define RMG_THREADS_H 1


/* Thread control structures */
typedef struct
{

    /* Thread ID number assigned by us */
    int tid;

    /* Assigned job */
    int job;

    /* With the complex option this lets the threads know which k-point is
     * currently being worked on in ortho and subdiag. */
    int kidx;

    /* Pointer to current state assigned to the thread when used in sections that process a single state */
    STATE *sp;

    rmg_double_t *vtot;

    // Pointers to special args
    void *p1;
    void *p2;
    int ion;        // Used for threaded beta_xpsi
    int nion;       // Used for threaded beta_xpsi
    rmg_double_t *sintR;    // Used for threaded beta_xpsi
    rmg_double_t *sintI;    // Used for threaded beta_xpsi
    rmg_double_t *weiptr;   // Used for threaded beta_xpsi
    int kpt;    // Used for threaded beta_xpsi
} SCF_THREAD_CONTROL;

/* Extern declarations for thread control structures */
extern SCF_THREAD_CONTROL thread_control[];

#endif
