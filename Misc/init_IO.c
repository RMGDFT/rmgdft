/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/***** RMG: Common/init_IO.c *****
 * NAME
 *   Ab initio real space multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *   Copyright (C) 2009  Frisco Rose, Jerzy Bernholc
 * FUNCTION
 *   void init_IO( int argc, char **argv )
 *   Initializes settings and creates directory structures for ouput logging
 *   Make each run image manage its own directory of input/output
 * INPUTS
 *   argc and argv from main
 * OUTPUT
 *   none
 * PARENTS
 *   main.c
 * CHILDREN
 *   init_pe.c read_pseudo.c
 * SOURCE
 */



#include "grid.h"
#include "const.h"
#include "params.h"
#include "rmgtypes.h"
#include "rmg_alloc.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "macros.h"

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>
#include "hybrid.h"

void init_IO (int argc, char **argv)
{

    int i, npes, worldpe, image, status, lognum = 0, provided=0, retval;
    char workdir[MAX_PATH], logname[MAX_PATH], basename[MAX_PATH], *quantity, *extension, *endptr;
    struct stat buffer;
    time_t timer;

    /* Set start of program time */
    timer = time (NULL);
    MPI_Init_thread(&argc, &argv, ct.mpi_threadlevel, &provided);

    /* get this cores mpi rank */
    MPI_Comm_rank (MPI_COMM_WORLD, &worldpe);
    pct.worldrank = worldpe;

    /* get total mpi core count */
    MPI_Comm_size (MPI_COMM_WORLD, &npes);
    pct.total_npes = npes;

    if(argc == 2)
    {
        read_init(argv[1]);
    }
    else
        read_init("ctrl_init.dat");


    init_pestr ();

    snprintf (ct.cfile, MAX_PATH, "%s%s", pct.image_path[pct.thisimg], pct.image_input[pct.thisimg]);
    snprintf (ct.basename, MAX_PATH, "%s%s", pct.image_path[pct.thisimg], pct.image_input[pct.thisimg]);


    read_control(ct.cfile);
    set_rank(pct.gridpe);

    init_HYBRID_MODEL(NPES, pct.gridpe, ct.THREADS_PER_NODE, pct.grid_comm);

    /* if logname exists, increment until unique filename found */
    if (pct.imgpe == 0)
    {
        snprintf (logname, MAX_PATH, "%s.log", ct.basename);

        int name_incr;
        name_incr = filename_increment(logname);
        snprintf (ct.logname, MAX_PATH, "%s.%02d", logname, name_incr);

        /* open and save logfile handle, printf is stdout before here */
        ct.logfile = fopen(ct.logname, "w");
    }
    else {
        ct.logfile = fopen("/dev/null", "w");
    }

    MPI_Comm_size (pct.img_comm, &status);
    printf ("\nRMG run started at GMT %s", asctime (gmtime (&timer)));
    printf ("\nRMG running with %d images and %d images per node.\n", pct.images, ct.images_per_node);
    printf ("\nRMG running in message passing mode with %d procs for this image.", status);

    /* Read in our pseudopotential information */
    ReadPseudo(ct.num_species, &ct.sp[0]);
   // read_pseudo();


#if GPU_ENABLED
    cudaDeviceReset();
    cudaSetDeviceFlags(cudaDeviceScheduleSpin);
    if( CUDA_SUCCESS != cuInit( 0 ) ) {
        fprintf(stderr, "CUDA: Not initialized\n" ); exit(-1);
    }
    if( CUDA_SUCCESS != cuDeviceGet( &ct.cu_dev, 0 ) ) {
        fprintf(stderr, "CUDA: Cannot get the device\n"); exit(-1);
    }
    cudaSetDevice(ct.cu_dev);
    if( CUBLAS_STATUS_SUCCESS != cublasInit( ) ) {
        fprintf(stderr, "CUBLAS: Not initialized\n"); exit(-1);
    }
    if( CUBLAS_STATUS_SUCCESS != cublasCreate(&ct.cublas_handle) ) {
        fprintf(stderr, "CUBLAS: Handle not created\n"); exit(-1);
    }
#if MAGMA_LIBS
magma_init();
//magmablasSetKernelStream(ct.cuda_stream);
#endif



#endif

    // This is placed down here since the IO is not setup yet when provided is obtained above.
    if(provided < ct.mpi_threadlevel) {

        printf("Thread support requested = %d but only %d provided. Terminating.\n", ct.mpi_threadlevel, provided);
        MPI_Finalize();
        exit(0);

    }
    printf("Running with thread level = %d\n", provided);
    fflush(NULL);

    // Allocate storage for trade_images and global sums routines
    init_TradeImages();
    set_MPI_comm(pct.grid_comm);

    init_global_sums();


}
