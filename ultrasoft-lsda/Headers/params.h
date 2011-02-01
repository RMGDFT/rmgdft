/************************** SVN Revision Information **************************
 **    $Id$    **
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


/* Makefile generated compile time parameters */
#include "arch.h"


/* Global grid sizes set in the Makefile */
#include "make_conf.h"


/*************************************************************
 *************************************************************

    All compile time parameters are now set in the top-level
    Makefile. Don't touch anything in this file unless you
    know what you are doing.

 *************************************************************
 *************************************************************/


/* This is the maximum number of non-local projectors for any species
 * just an upper limit on dynamically allocated data structures. */
#define         MAX_NL  	18
#define         MAX_NB  	6


/* Maximum number of species -- just an upper limit on
 * dynamically allocated data structures. */
#define         MAX_SPECIES     64


/* Maximum number of ions -- just an upper limit on
 * dynamically allocated data structures. */
#define         MAX_IONS        512


/* Maximum l-value (Angular momentum channel for pseudotentials) */
#define         MAX_L       	4


// Default character buffer size
#define		MAX_CHAR	255

/* Maximum pathname length for input file names */
#define         MAX_PATH        MAX_CHAR

/* Maximum number of points in pseudopotential radial grid */
#define         MAX_RGRID	(1000)



/* Size of the linear interpolation grid for the potentials */
#define         MAX_LOCAL_LIG	(9000)


/* Size of the linear interpolation grid for the potentials */
#define         MAX_QLIG       	(9000)


/* Dimensions of shmem workspace */
#define		MAX_PWRK	(50000)


/* The number of possible point symmetries */
#define MAX_SYMMETRY	48

#ifdef SMP
#define 	MAX_THREADS	1000
#else
#define 	MAX_THREADS	1
#endif





#ifdef SMP

#  define  	PE_X  		1
#  define  	PE_Y  		1
#  define  	PE_Z  		1

#elif MPI

#  ifndef PE_X
    /* Processor topology */
#    if   NPES == 1
#      define		PE_X		1
#      define		PE_Y	        1
#      define		PE_Z	        1
#    elif NPES == 2
#      define		PE_X		2
#      define		PE_Y	        1
#      define		PE_Z	        1
#    elif NPES == 4
#      define		PE_X		2
#      define		PE_Y	        2
#      define		PE_Z	        1
#    elif NPES == 8
#      define		PE_X            2
#      define		PE_Y            2
#      define		PE_Z            2
#    elif NPES == 16
#      define		PE_X            2
#      define		PE_Y            4
#      define		PE_Z            2
#    elif NPES == 32
#      define		PE_X            4
#      define		PE_Y            4
#      define		PE_Z            2
#    elif NPES == 48
#      define		PE_X            4
#      define		PE_Y            2
#      define		PE_Z            6
#    elif NPES == 64
#      define		PE_X            4
#      define		PE_Y            4
#      define		PE_Z            4
#    elif NPES == 128
#      define		PE_X            8
#      define		PE_Y            4
#      define		PE_Z            4
#    elif NPES == 256
#      define		PE_X            8
#      define		PE_Y            8
#      define		PE_Z            4
#    elif NPES == 512
#      define		PE_X            8
#      define		PE_Y            8
#      define		PE_Z            8
#    endif
#  endif

#endif


#ifdef SMP
#  define         PHYS_PES        1
#else
#  define         PHYS_PES        (NPES)
#endif






/* Fine grid size on each processor */
#define         FNX_GRID        (NX_GRID * FG_NX)
#define         FNY_GRID        (NY_GRID * FG_NY)
#define         FNZ_GRID        (NZ_GRID * FG_NZ)

#define         PX0_GRID        (NX_GRID / PE_X)
#define         PY0_GRID        (NY_GRID / PE_Y)
#define         PZ0_GRID        (NZ_GRID / PE_Z)

#define         P0_BASIS        (PX0_GRID * PY0_GRID * PZ0_GRID)

#define         FPX0_GRID       (PX0_GRID * FG_NX)
#define         FPY0_GRID       (PY0_GRID * FG_NY)
#define         FPZ0_GRID       (PZ0_GRID * FG_NZ)

#define         FP0_BASIS        (FPX0_GRID * FPY0_GRID * FPZ0_GRID)




/******/
