#ifndef RMG_THREADS_H
#define RMG_THREADS_H 1


#ifdef USE_NUMA
    #include <numa.h>
#endif


/* Thread control structures */
typedef struct
{

    /* Assigned job */
    int job;

    // basetag which is transferred by the project specific thread routine
    // into the BaseThreadControl structure.
    int basetag;

    // vcycle during multigrid iterations
    int vcycle;

    /* Pointer to current state assigned to the thread when used in sections that process a single state */
    void *sp;

    double *vtot;

    double *rho_neutral;
    int boundaryflag;

    // Pointers to special args
    void *p1;
    void *p2;
    void *p3;
    void *p4;
    void *nv;         // Non-local operator applied to a specific wavefunction
    void *ns;         // S-operator applied to a specific wavefunction
    void *Bns;        // Bapplied to S-operator applied to a specific wavefunction
    double eig;            // Used for Davidson preconditioner
    double avg_potential;  // Used for Davidson preconditioner
    double fd_diag;        // Used for Davidson preconditioner
    int istate;
} SCF_THREAD_CONTROL;

// Called from main to setup thread tasks
void QueueThreadTask(int tid, SCF_THREAD_CONTROL &task);
// Called from threads to get what they are supposed to do
bool PopThreadTask(int tid, SCF_THREAD_CONTROL &task);
// Called from main to terminate all threads
void RmgTerminateThreads(void);
#endif
