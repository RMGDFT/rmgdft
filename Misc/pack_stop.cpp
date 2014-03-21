/************************** SVN Revision Information **************************
 **    $Id: pack_stop.c 2093 2014-02-06 01:02:40Z ebriggs $    **
******************************************************************************/

/****f* QMD-MGDFT/pack_stop.c *****
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
 *   void pack_stop(rmg_double_t *sg, rmg_double_t *pg, int dimx, int dimy, int dimz)
 *   pack smoothing grids sg to processor grids pg
 * INPUTS
 *   sg[(dimx+2)*(dimy+2)*(dimz+2)]: array with images
 *   dimx, dimy, dimz: dimensions of array 
 * OUTPUT
 *   pg[dimx*dimy*dimz]: processor grid, its value are copied from sg
 * PARENTS
 *   get_vh.c
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include "make_conf.h"
#include "rmgtypes.h"
#include "auxiliary.h"
#include "BlasWrappers.h"

using namespace std;

template void CPP_pack_stop<double>(double*, double*, int, int, int);
template void CPP_pack_stop<float>(float*, float*, int, int, int);

template <typename RmgType>
void CPP_pack_stop (RmgType * sg, RmgType * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
    int incx, incy, incxs, incys;
    int ione = 1;
    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);


    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            QMD_copy (dimz, &sg[ixh * incxs + iyh * incys + 1], ione, &pg[ix * incx + iy * incy],
                       ione);

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_stop */


extern "C" void pack_stop (rmg_double_t * sg, rmg_double_t * pg, int dimx, int dimy, int dimz)
{
    CPP_pack_stop<rmg_double_t> (sg, pg, dimx, dimy, dimz);
}

extern "C" void pack_stop_f (rmg_float_t * sg, rmg_float_t * pg, int dimx, int dimy, int dimz)
{
    CPP_pack_stop<rmg_float_t> (sg, pg, dimx, dimy, dimz);
}
