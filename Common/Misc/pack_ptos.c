/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/pack_ptos.c *****
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
 *   void pack_ptos(REAL *sg, REAL *pg, int dimx, int dimy, int dimz)
 *   pack processor grids pg to smoothing grids sg
 * INPUTS
 *   pg[dimx*dimy*dimz]: original array 
 *   dimx, dimy, dimz: dimensions of the array
 * OUTPUT
 *   sg[(dimx+2)*(dimy+2)*(dimz+2)]: its elements except for images are copied from pg
 * PARENTS
 *   get_ke.c get_vh.c init_wf.c mg_eig_state.c subdiag_mpi.c subdiag_smp.c xcgga.c
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include "md.h"
#include <float.h>
#include <math.h>



void pack_ptos(REAL * sg, REAL * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
    int incx, incy, incxs, incys;
    int ione = 1;
    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);

    for (ix = 0; ix < (dimx + 2) * (dimy + 2) * (dimz + 2); ix++)
        sg[ix] = 0.0;

    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            QMD_scopy(dimz, &pg[ix * incx + iy * incy], ione,
                      &sg[ixh * incxs + iyh * incys + 1], ione);

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_ptos */

/******/
