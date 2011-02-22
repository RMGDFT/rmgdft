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
 *   REAL real_min_all(REAL x)
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



REAL real_min_all (REAL x)
{

#ifndef MPI

    return x;

#else
    REAL inreg;
    REAL outreg;

    inreg = x;

    MPI_Allreduce (&inreg, &outreg, 1, MPI_DOUBLE, MPI_MIN, pct.grid_comm);

    return outreg;


#endif


}                               /* end real_sum_all */


REAL real_min_all_spin (REAL x)
{

#ifndef MPI

    return x;

#else
    REAL inreg;
    REAL outreg;

    inreg = x;

    MPI_Allreduce (&inreg, &outreg, 1, MPI_DOUBLE, MPI_MIN, pct.img_comm);

    return outreg;


#endif


}                               /* end real_sum_all */

/******/
