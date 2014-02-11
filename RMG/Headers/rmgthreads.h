#ifndef RMG_THREADS_H
#define RMG_THREADS_H 1


/* Thread control structures */
typedef struct
{

    /* Thread ID number assigned by us */
    int tid;

#if GPU_ENABLED
    // Cuda device context
    cudaStream_t cstream;
    rmg_double_t *gpu_host_temp1;
    rmg_double_t *gpu_host_temp2;
#endif

    /* Thread identifier from pthread_self. Needed to send signals */
    pthread_t pthread_tid;

    /* Assigned job */
    int job;

    /* Synchronization semaphore */
    sem_t sync;

    /* These volatiles are used as synchronization variables for the threads */
    volatile int start;

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
