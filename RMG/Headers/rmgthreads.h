#ifndef RMG_THREADS_H
#define RMG_THREADS_H 1


/* Thread control structures */
typedef struct
{

    /* Thread ID number assigned by us */
    int tid;

    /* MPI communicator for use by this thread */
    MPI_Comm grid_comm;

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

    /* Pointer to state array used by each thread */
    STATE *my_states;

    /* Local variable -- summed to obtain total charge for all orbitals */
    rmg_double_t tcharge;

    /* Spacial offset for the thread */
    int offset;

    /* Points to base of distributed storage array for this thread */
    rmg_double_t *base_mem;

    /* Points to base of distributed scratch array for this thread */
    rmg_double_t *scratch1;

    /* Number of points per wavefunction in the distributed storage array */
    int numpt;

    /* leading dimension of the distributed wave function storage array */
    int lda;

    /* Local copies of eigenvalues and occupations for this thread */
    rmg_double_t *eigs;
    rmg_double_t *occs;

    /* Force contributions computed by this thread */
    rmg_double_t force[MAX_IONS][3];

    /* Pointer to dynamically allocated arrays of size ct.num_states*ct.num_states */
    /* that is required in ortho. Each thread has it's own copy */
    rmg_double_t *darr;
    rmg_double_t *barr;


    /* The same array as referenced by darr but this copy is 
     *allocated in the main program rather than in one of the threads.
     */
    rmg_double_t *farr;


    rmg_double_t *rho;
    rmg_double_t *rhocore;
    rmg_double_t *vtot;
    rmg_double_t *vnuc;

    /* Pointers to the non-local potential index list 
     *and to the projectors themselves */
    int *nlindex;
    rmg_double_t *projectors;

    // Pointers to special args
    void *p1;
    void *p2;
    rmg_double_t *trade_buf;// Used by trade_images
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
