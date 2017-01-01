#ifndef RMG_THREADS_H
#define RMG_THREADS_H 1


#ifdef USE_NUMA
    #include <numa.h>
#endif

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

    // vcycle during multigrid iterations
    int vcycle;

    /* Pointer to current state assigned to the thread when used in sections that process a single state */
    void *sp;

    double *vtot;

    double *rho_neutral;
    double rms_target;
    int min_sweeps, max_sweeps,mg_level;
    int boundaryflag;

    // Pointers to special args
    void *p1;
    void *p2;
    void *p3;
    void *p4;
    int ion;          // Used for threaded beta_xpsi
    int nion;         // Used for threaded beta_xpsi
    double *sintR;    // Used for threaded beta_xpsi
    double *sintI;    // Used for threaded beta_xpsi
    double *weiptr;   // Used for threaded beta_xpsi
    void *nv;         // Non-local operator applied to a specific wavefunction
    void *ns;         // S-operator applied to a specific wavefunction
    void *Bns;        // Bapplied to S-operator applied to a specific wavefunction
    int kpt;          // Used for threaded beta_xpsi
    double eig;            // Used for Davidson preconditioner
    double avg_potential;  // Used for Davidson preconditioner
    double fd_diag;        // Used for Davidson preconditioner
    int istate;
} SCF_THREAD_CONTROL;

/* Extern declarations for thread control structures */
extern SCF_THREAD_CONTROL thread_control[];

#endif
