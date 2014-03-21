/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/real_sum_all.c *****
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
 *   rmg_double_t real_min_all(rmg_double_t x)
 *   Finds minimum over all processors.
 * INPUTS
 *   x: defined in each processor
 * OUTPUT
 *   mimimum over all processors is returned to each processor
 * PARENTS
 * ??
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include "main.h"



rmg_double_t real_min_all (rmg_double_t x, MPI_Comm comm)
{

    rmg_double_t inreg;
    rmg_double_t outreg;

    inreg = x;

    MPI_Allreduce (&inreg, &outreg, 1, MPI_DOUBLE, MPI_MIN, comm);

    return outreg;

}                               /* end real_sum_all */

