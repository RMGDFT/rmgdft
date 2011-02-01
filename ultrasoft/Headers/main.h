/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*** QMD-MGDFT/main.h *****
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

                        main.h

    Header file for molecular dynamics code.



*/

#define REAL_SPACE 1


/* stdio for file handle argument type */
#include <stdio.h>
#include <stdlib.h>

/* Use predefined booleans */
#include <stdbool.h>

/* Size of floating point variables used in QMD */
#define     REAL    double

/* Version information */
#include    "version.h"

/* Compile time parameters */
#include    "params.h"

/* Constants and symbolic definitions */
#include    "const.h"

/* include the mpi wrapper */
#include    "my_mpi.h"

/* include scalapack wrapper */
#include "my_scalapack.h"

/* fourier transformation structure definitions*/
#include    "fftw.h"


/* Header file for blas routines */
#include "blas.h"


/* this declares macros and external functions for memory allocation */
#include "salloc.h"


/* other general macros */
#include "macros.h"

/* routines for input parsing */
#include "input.h"

/* Custom types used in the code*/
#include "typedefs.h"

/*Prototypes for function calls*/
#include "prototypes.h"

/******/
