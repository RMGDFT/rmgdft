/****f* QMD-MGDFT/params.h *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   
 * INPUTS
 *
 * OUTPUT
 *  
 * PARENTS
 *
 * CHILDREN
 * 
 * SEE ALSO
 *  
 * SOURCE
 */

/*

                            params.h


    Compile time parameters for md program.


*/

#ifndef RMG_params_H
#define RMG_params_H 1





/*************************************************************
 *************************************************************

    All compile time parameters are now set in the top-level
    Makefile. Don't touch anything in this file unless you
    know what you are doing.

 *************************************************************
 *************************************************************/

/* Maximum l-value (Angular momentum channel for pseudotentials) */
#define         MAX_L       	4

// Default character buffer size
#define		MAX_CHAR	255

/* Maximum pathname length for input file names */
#define         MAX_PATH        MAX_CHAR

/* Maximum image count */
#define         MAX_IMGS        99

/* Maximum number of points in pseudopotential radial grid */
#define         MAX_RGRID	(2000)

/* Maximum number of points to use in the log interpolation grid */
#define         MAX_LOGGRID      (2000)

/* Starting radius for the log interpolation grid */
#define         LOGGRID_START    (0.00001)

/* Log interpolation grid mesh parameter */
#define         LOGGRID_MESHPARM  (1.0075)


/* Dimensions of shmem workspace */
#define		MAX_PWRK	(1024)


/* The number of possible point symmetries */
#define MAX_SYMMETRY	48



/* Size of the linear interpolation grid for the qfunction */
/* Max. number of states localized on the same ion */
#define MAX_LOC_ST      20

#define MAX_BLOCKS      20

/* Restart files are compressed using zfp compressor and this accuracy factor */
#define RESTART_TOLERANCE  1.0e-09

#if CUDA_ENABLED || HIP_ENABLED

    #define MAX_GPU_DEVICES 16

#endif

#endif
