/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/app6_del2.c *****
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
 *   void app6_del2(rmg_double_t *f, rmg_double_t *work, int dimx, int dimy, int dimz)
 *   Apply the sixth order Laplacian operator on a orthorhombic grid
 * INPUTS
 *   S0_GRID *f:  real array of (dimx+2) * (dimy+2) * ((dimz+2)
 *   see main.h for defination of unions S0_GRID and P0_GRID
 * OUTPUT
 *   P0_GRID *work: real array of dimx * dimy * dimz
 * PARENTS
 *   get_ke.c xcgga.c
 * CHILDREN
 *   trade_image2.c
 * SOURCE
 */



#include "main.h"
#include "prototypes_on.h"
#include <float.h>
#include <math.h>


void app6_del2(rmg_double_t * f, rmg_double_t * work, int dimx, int dimy, int dimz,
               rmg_double_t hxgrid, rmg_double_t hygrid, rmg_double_t hzgrid)
{

    int iz, ix, iy;
    rmg_double_t h2, t0, t1x, t2x;
    rmg_double_t t1y, t2y;
    rmg_double_t t1z, t2z;
    double t3x,t3y,t3z;
    rmg_double_t *dum2;
    int ixs, iys, ix1, iy1;

#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc();
#endif

    ixs = (dimy + 6) * (dimz + 6);
    iys = (dimz + 6);
    ix1 = dimy * dimz;
    iy1 = dimz;

    my_malloc_init( dum2, (dimx + 6) * (dimy + 6) * (dimz + 6), rmg_double_t );
    trade_imagesx(f, dum2, dimx, dimy, dimz, 3, CENTRAL_FD);

    h2 = hxgrid * hxgrid * get_xside() * get_xside();
    t0 = -49.0 / (18.0 * h2);
    t1x =  3.0 / ( 2.0 * h2);
    t2x = -3.0 / (20.0 * h2);
    t3x =  1.0 / (90.0 * h2);

    h2 = hygrid * hygrid * get_yside() * get_yside();
    t0 -= 49.0 / (18.0 * h2);
    t1y =  3.0 / ( 2.0 * h2);
    t2y = -3.0 / (20.0 * h2);
    t3y =  1.0 / (90.0 * h2);

    h2 = hzgrid * hzgrid * get_zside() * get_zside();
    t0 -= 49.0 / (18.0 * h2);
    t1z =  3.0 / ( 2.0 * h2);
    t2z = -3.0 / (20.0 * h2);
    t3z =  1.0 / (90.0 * h2);



    for (ix = 3; ix < dimx + 3; ix++)
    {

        for (iy = 3; iy < dimy + 3; iy++)
        {

            for (iz = 3; iz < dimz + 3; iz++)
            {

                work[(ix - 3) * ix1 + (iy - 3) * iy1 + iz - 3] =
                    t0 * dum2[ix * ixs + iy * iys + iz] +
                    t1x * dum2[(ix - 1) * ixs + iy * iys + iz] +
                    t1x * dum2[(ix + 1) * ixs + iy * iys + iz] +
                    t2x * dum2[(ix - 2) * ixs + iy * iys + iz] +
                    t2x * dum2[(ix + 2) * ixs + iy * iys + iz] +
                    t3x * dum2[(ix - 3) * ixs + iy * iys + iz] +
                    t3x * dum2[(ix + 3) * ixs + iy * iys + iz] +
                    t1y * dum2[ix * ixs + (iy - 1) * iys + iz] +
                    t1y * dum2[ix * ixs + (iy + 1) * iys + iz] +
                    t2y * dum2[ix * ixs + (iy - 2) * iys + iz] +
                    t2y * dum2[ix * ixs + (iy + 2) * iys + iz] +
                    t3y * dum2[ix * ixs + (iy - 3) * iys + iz] +
                    t3y * dum2[ix * ixs + (iy + 3) * iys + iz] +
                    t1z * dum2[ix * ixs + iy * iys + iz - 1] +
                    t1z * dum2[ix * ixs + iy * iys + iz + 1] +
                    t2z * dum2[ix * ixs + iy * iys + iz - 2] +
                    t2z * dum2[ix * ixs + iy * iys + iz + 2] +
                    t3z * dum2[ix * ixs + iy * iys + iz - 3] +
                    t3z * dum2[ix * ixs + iy * iys + iz + 3];

            }                   /* end for */
        }                       /* end for */
    }                           /* end for */

#if MD_TIMERS
    time2 = my_crtc();
    rmg_timings(APPCIL_TIME, (time2 - time1));
#endif

    my_free(dum2);

}                               /* end app6_del2 */


/******/
