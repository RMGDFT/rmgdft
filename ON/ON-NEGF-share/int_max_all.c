/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/int_max_all.c *****
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
 *   int int_max_all(int x)
 *   Performs a scalar sum over all processors.
 * INPUTS
 *   x: defined in each processor
 * OUTPUT
 *   sum over all processors is returned to each processor
 * PARENTS
 *   get_ke.c get_rho.c get_te.c get_vh.c getpoi_bc.c init_nuc.c
 *   lforce.c mg_eig_state.c norm_psi.c 
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include "main.h"
#include "prototypes_on.h"
#include <float.h>
#include <math.h>


int int_max_all(int x)
{

    int inreg;
    int outreg;

    inreg = x;
#if SERIAL

    outreg = x;

#else

    MPI_Allreduce(&inreg, &outreg, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);

#endif

    return outreg;


}                               /* end real_max_all */


/******/
