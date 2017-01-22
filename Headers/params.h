/************************** SVN Revision Information **************************
 **    $Id: params.h 1870 2013-02-05 13:58:02Z ebriggs $    **
******************************************************************************/
 
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


#include "make_conf.h"


/*************************************************************
 *************************************************************

    All compile time parameters are now set in the top-level
    Makefile. Don't touch anything in this file unless you
    know what you are doing.

 *************************************************************
 *************************************************************/

/* Maximum number of PE's */
#define         MAX_PES        32768

/* This is the maximum number of non-local projectors for any species
 * just an upper limit on dynamically allocated data structures. */
#define         MAX_NL  	38
#define         MAX_NB  	9


/* Maximum number of species -- just an upper limit on
 * dynamically allocated data structures. */
#define         MAX_SPECIES     64


/* Maximum number of ions -- just an upper limit on
 * dynamically allocated data structures. */
#define         MAX_IONS        5120

/* Maximum number of atoms whose projectors have overlap with current processor */
#define MAX_NONLOC_IONS         1024

/* Maximum number of processors with which the current processor will communicate to calculate <beta|psi> */
#define MAX_NONLOC_PROCS        250

/* Maximum l-value (Angular momentum channel for pseudotentials) */
#define         MAX_L       	4

// Default character buffer size
#define		MAX_CHAR	255

/* Maximum pathname length for input file names */
#if (defined(_WIN32) || defined(_WIN64))
    #include <windows.h>
#else
    #define         MAX_PATH        MAX_CHAR
#endif

/* Maximum image count */
#define         MAX_IMGS        99

/* Maximum cpu grids, now spin enabled -> 2.
 * This will need to be dynamic if we parallelize over k-points */
#define         MAX_GRIDS        2

/* Maximum number of points in pseudopotential radial grid */
#define         MAX_RGRID	(2000)



/* Size of the linear interpolation grid for the potentials */
#define         MAX_LOCAL_LIG	(9000)

/* Maximum number of points to use in the log interpolation grid */
#define         MAX_LOGGRID      (2000)

/* Starting radius for the log interpolation grid */
#define         LOGGRID_START    (0.00001)

/* Log interpolation grid mesh parameter */
#define         LOGGRID_MESHPARM  (1.0075)


/* Size of the linear interpolation grid for the potentials */
#define         MAX_QLIG       	(9000)

/* Dimensions of shmem workspace */
#define		MAX_PWRK	(1024)


/* The number of possible point symmetries */
#define MAX_SYMMETRY	48

#define         MAX_OPTS  16



/* Number of k-points -- just an upper limit on
 * dynamically allocated data structures. */
#define         MAX_KPTS           32


/* Size of the linear interpolation grid for the qfunction */
/* Max. number of states localized on the same ion */
#define MAX_LOC_ST      20

#define MAX_BLOCKS      20


#if GPU_ENABLED

    #define MAX_GPU_DEVICES 16

#endif

#endif
