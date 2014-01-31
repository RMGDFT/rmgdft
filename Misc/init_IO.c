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
#include "main.h"
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>
#if HYBRID_MODEL
#include <pthread.h>
#include "hybrid.h"
#endif

#if PAPI_PERFMON
#undef kill
#include <papi.h>
PAPI_option_t papi_opts;
#endif

void init_IO (int argc, char **argv)
{

    int i, npes, worldpe, image, status, lognum = 0, provided=0, retval;
    char workdir[MAX_PATH], logname[MAX_PATH], basename[MAX_PATH], *quantity, *extension, *endptr;
    struct stat buffer;
    time_t timer;

#if PAPI_PERFMON
    PAPI_option_t options;
#endif

    /* Set start of program time */
    timer = time (NULL);
#if HYBRID_MODEL
    MPI_Init_thread(&argc, &argv, ct.mpi_threadlevel, &provided);
#else

    /* Initialize MPI, we need it for error_handler, amongst others */
    MPI_Init (&argc, &argv);

#endif

    /* get this cores mpi rank */
    MPI_Comm_rank (MPI_COMM_WORLD, &worldpe);

    /* get total mpi core count */
    MPI_Comm_size (MPI_COMM_WORLD, &npes);

    /* Define a default output stream, gets redefined to log file later */
    ct.logfile = stdout;

    /* 1st argument must exist as a file, use it as the control file or die */
    if ((status = stat (argv[1], &buffer)) == -1)
        error_handler ("Command line argument \"%s\" is malformed, file D.N.E.", argv[1]);
    strncpy (ct.cfile, argv[1], MAX_PATH);

    /* test to see if argument 1 determines a multi-image run, fail gracefully if not */
    pct.images = 1;
    strncpy (basename, argv[1], MAX_PATH);
    if ((extension = rindex (basename, '.')) != NULL)
        *extension++ = '\0';
    if ((quantity = rindex (basename, '.')) != NULL)
    {
        *quantity++ = '\0';
        pct.images = (int) strtol (quantity, &endptr, 10);

        /*Check if this is *just* a number, not some numerical naming convention*/
        if (quantity == endptr)
            pct.images = 1;
    }

    dprintf ("\nRunning with %d images", pct.images);
    if ( pct.images > MAX_IMGS )
        error_handler ("Multi-image input file %s asks for more images (%d) than MAX_IMGS in params.h.", argv[1], pct.images);

    /* Setup image number that this core belongs to */
    image = npes / pct.images;
    NPES = image;   // NPES is the per image number of PES

    if (image * pct.images != npes)
        error_handler ("Total MPI processes (%d) must be a multiple of the number of images (%d) in this run.", npes, pct.images);

    // Only works for 2 images stacked vertically.
    if(ct.images_per_node > 1) {
        int image_x, image_xs;
        image_xs = ct.images_per_node;
        
        if(worldpe % ct.images_per_node) {
            // image_y=1
            image_x = (worldpe - 1) / (2 * NPES);
            pct.thisimg = image_xs * image_x + 1;
        }
        else {
            // image_y=0
            image_x = worldpe / (2 * NPES);
            pct.thisimg = image_xs * image_x;
        }
    }
    else {
        pct.thisimg = worldpe / image;
    }

    if (pct.images > 1)
        snprintf (ct.cfile, MAX_PATH, "%s.%02d.rmg", basename, pct.thisimg + 1);

    //error_handler("message a"); 
    /* PE(MPI) initialization, need mpi groups defined before logfile is initialized */
    init_pe ( image );

    //error_handler("message b");  

    if (pct.images > 1)
    {
        if (strcmp (extension, "rmg") != 0)
            error_handler 
		    ("Multi-image input file %s does not end with proper \"rmg\" extension.", argv[1]);

        /* logfile name is based on input file and this images group number */
        /* if second command line argument exists, use it as a basename for output */
        if (argc == 3)
        {
            snprintf (ct.basename, MAX_PATH, "%s.%02d", argv[2], pct.thisimg + 1);
        }
        else
        {
            snprintf (ct.basename, MAX_PATH, "%s.%02d", basename, pct.thisimg + 1);
        }

        /* every image has it own output/working directory */
        sprintf (workdir, "image.%02d", pct.thisimg + 1);
        if (pct.imgpe == 0)
        {
            if (status = stat ( workdir, &buffer ) == 0)
            {
                if ( !S_ISDIR( buffer.st_mode))
                    error_handler 
			    ("Found %s, that is not a directory as required for multi-instance run!", workdir);
            }
            else
            {
                if ((status = mkdir (workdir, S_IRWXU)) == -1)
                    error_handler
                        ("Unable to create image directory \"%s\", mkdir returned %d, check file system for permissions/diskfull.",
                         workdir, status);
            }
        }
        /* Make sure that process 0 has successfully verified/created the image directory */
        my_barrier ();

        if ((status = chdir (workdir)) == -1)
            error_handler
                ("Unable to cd to image directory \"%s\", check file system for permissions/diskfull.",
                 workdir);
    }
    else
    {
        /* Not multi-image run, so set output file accordingly */
        /* if second command line argument exists, use it as a basename for output */
        if (argc == 3)
        {
            snprintf (ct.basename, MAX_PATH, "%s", argv[2]);
        }
        else
        {
            /* Reset basename to be input file name (except rmg extension) */
            strncpy (basename, argv[1], MAX_PATH);
            if ((extension = rindex (basename, '.')) != NULL)
            {
                    *extension++ = '\0';
                if (strcmp (extension, "rmg") == 0)
                    snprintf (ct.basename, MAX_PATH, "%s", basename);
                else
                    snprintf (ct.basename, MAX_PATH, "%s.%s", basename, extension);
            }
            else
            {
                snprintf (ct.basename, MAX_PATH, "%s", basename);
            }
        }
    }

    /* if logname exists, increment until unique filename found */
    if (pct.imgpe == 0)
    {
	snprintf (logname, MAX_PATH, "%s.log", ct.basename);
        while ((status = stat (logname, &buffer)) != -1)
        {
            strncpy (basename, logname, MAX_PATH);
            if (++lognum > 99)
                error_handler
			("You have over 100 logfiles, you need to think of a better job scenario!\n");
            if ((extension = rindex (basename, '.')) != NULL)
                *extension = '\0';
            if ((quantity = rindex (basename, '.')) != NULL)
                if ( quantity[1] == 'r' )
                {
                    /* throw away image number, just satisfies function return */
                     image = (int) strtol (&quantity[2], &endptr, 10);

                    /*Check if this is a log rename number, additionally it should only be 2 digits long*/
                    if (extension == endptr)
                        *quantity = '\0';
                }
            snprintf (ct.basename, MAX_PATH, "%s.r%02d", basename, lognum );
            snprintf (logname, MAX_PATH, "%s.log", ct.basename);
        }

        /* open and save logfile handle, printf is stdout before here */
        my_fopen (ct.logfile, logname, "w");
    }

    MPI_Comm_size (pct.img_comm, &status);
    printf ("\nRMG run started at GMT %s", asctime (gmtime (&timer)));
    printf ("\nRMG running in message passing mode with %d procs.", status);
    if (pct.images == 1)
        printf ("\nRMG running in a single-image mode\n");
    else
        printf ("\nRMG running in a multi-image mode with %d images and %d images per node.\n", pct.images, ct.images_per_node);

    /* Read in our pseudopotential information */
    read_pseudo ();

    // If papi performance monitoring is desired initialize the library
#if PAPI_PERFMON
    if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
        error_handler("Problem with papi support. Exiting.\n");
    }
    if (PAPI_thread_init(pthread_self) != PAPI_OK) {
        error_handler("Problem with papi thread support. Exiting.\n");
    }

#if 1
    ct.EventSet = PAPI_NULL;
    if (PAPI_create_eventset(&ct.EventSet) != PAPI_OK) {
         error_handler ("Cannot create PAPI event set.\n");
    }
    if (PAPI_add_event(ct.EventSet, PAPI_FP_OPS) != PAPI_OK) {
         error_handler ("Cannot add PAPI event.\n");
    }

    PAPI_start(ct.EventSet);
     if(PAPI_set_granularity(PAPI_GRN_PROC) != PAPI_OK) {
         error_handler ("Cannot add PAPI granularity.\n");
     }

#endif

#pragma omp parallel for
    for(i = 0;i < THREADS_PER_NODE;i++) 
        Papi_init_omp_threads(i);
#endif

    // Allocate storage for trade_images and global sums routines
    trade_images (NULL, 0, 0, 0, NULL, 0);
    trade_images_f (NULL, 0, 0, 0, NULL, 0);
    trade_imagesx (NULL, NULL, 0, 0, 0, 0, 0);
    init_global_sums();

#if ASYNC_TRADES
  // set up memory sections for async trade images
  init_trade_imagesx_async();
#endif

#if GPU_ENABLED
  cudaDeviceReset();
  if( CUDA_SUCCESS != cuInit( 0 ) ) {
      fprintf(stderr, "CUDA: Not initialized\n" ); exit(-1);
  }
  if( CUDA_SUCCESS != cuDeviceGet( &ct.cu_dev, 0 ) ) {
      fprintf(stderr, "CUDA: Cannot get the device\n"); exit(-1);
  }
  cudaSetDevice(ct.cu_dev);
//  if( CUDA_SUCCESS != cuCtxCreate( &ct.cu_context, CU_CTX_SCHED_YIELD, ct.cu_dev ) ) {
//      fprintf(stderr, "CUDA: Cannot create the context\n"); exit(-1);
//  }
  cudaSetDevice(ct.cu_dev);
  if( CUBLAS_STATUS_SUCCESS != cublasInit( ) ) {
      fprintf(stderr, "CUBLAS: Not initialized\n"); exit(-1);
  }

#endif

#if HYBRID_MODEL
  if(provided < ct.mpi_threadlevel) {

      printf("Thread support requested = %d but only %d provided. Terminating.\n", ct.mpi_threadlevel, provided);
      MPI_Finalize();
      exit(0);

  }
  printf("Running with thread level = %d\n", provided);
  fflush(NULL);
  init_HYBRID_MODEL();
#endif
    return;
}
