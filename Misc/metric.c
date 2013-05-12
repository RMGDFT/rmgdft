/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/metric.c *****
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
 *   rmg_double_t metric(REAL *crystal)
 *   Calculate length (or distance) of a vector crystal(3)
 * INPUTS
 *   crystal: real array dimensioned 3
 * OUTPUT
 *   distance is returned
 * PARENTS
 *   cholesky.c get_nlop_d.c get_nlop_p.c get_nlop_s.c getpoi_bc.c init_nuc.c 
 *   init_wflcao.c lforce.c minimage.c moldyn.c nlforce_d.c nlforce_p.c nlforce_s.c
 *   subdiag_mpi.c subdiag_smp.c
 * CHILDREN
 *   to_crystal.c to_cartisian.c
 * SOURCE
 */
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

rmg_double_t metric (REAL * crystal)
{
    rmg_double_t cartesian[3];          /* cartesian coordinates of point */
    rmg_double_t distance;
    int ir;
    to_cartesian (crystal, cartesian);

    distance = 0.0;

    for (ir = 0; ir < 3; ir++)
        distance += cartesian[ir] * cartesian[ir];

    distance = sqrt (distance);

    return (distance);

}                               /* end metric */

/******/
