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


/* Typedefs */

#define REAL_SPACE 1


/* stdio for file handle argument type */
#include <stdio.h>
#include <stdlib.h>

/* Use predefined booleans */
#include <stdbool.h>

/* Typedefs */
#include    "rmgtypedefs.h"

/* Compile time parameters */
#include    "params.h"

/* include the mpi wrapper */
#include    "my_mpi.h"

/* include scalapack wrapper */
#include "my_scalapack.h"

/* Header file for blas routines */
#include "blas.h"


/* this declares macros and external functions for memory allocation */
#include "rmg_alloc.h"


/* other general macros */
#include "macros.h"

/* routines for input parsing */
#include "input.h"


/* Custom types used in the code*/
#include "typedefs.h"


#include "lbfgs.h"


/* Prototypes for function calls*/
#include "common_prototypes.h" 
#include "common_prototypes1.h" 

/* Some definitions needed for using compile time values for global grid sizes */
#include "fixed_dims.h"

#include "xc.h"

/* Include the  library of exchange and correlation (namely, Libxc) */
/* #include "../../Common/libxc/include/xc.h"   */


/******/
