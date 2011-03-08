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
 *   void app6_del2(REAL *f, REAL *work, int dimx, int dimy, int dimz)
 *   Apply the sixth order Laplacian operator on a orthorhombic grid
 * INPUTS
 *   S0_GRID *f:  real array of (dimx+2) * (dimy+2) * ((dimz+2)
 *   see md.h for defination of unions S0_GRID and P0_GRID
 * OUTPUT
 *   P0_GRID *work: real array of dimx * dimy * dimz
 * PARENTS
 *   get_ke.c xcgga.c
 * CHILDREN
 *   trade_image2.c
 * SOURCE
 */



#include "md.h"
#include <float.h>
#include <math.h>


void app4_del2(REAL * f, REAL * work, int dimx, int dimy, int dimz,
               REAL hxgrid, REAL hygrid, REAL hzgrid)
{

    int iz, ix, iy;
    REAL h2, t0, t1x, t2x;
    REAL t1y, t2y;
    REAL t1z, t2z;
    REAL *dum2;
    int ixs, iys, ix1, iy1;

#if MD_TIMERS
    REAL time1, time2;
    time1 = my_crtc();
#endif

    ixs = (dimy + 4) * (dimz + 4);
    iys = (dimz + 4);
    ix1 = dimy * dimz;
    iy1 = dimz;

    my_malloc_init( dum2, (dimx + 4) * (dimy + 4) * (dimz + 4), REAL );
    trade_images2(f, dum2, dimx, dimy, dimz);

    h2 = hxgrid * hxgrid * ct.xside * ct.xside;
    t0 = -30.0 / (12.0 * h2);
    t1x = 16.0 / (12.0 * h2);
    t2x = -1.0 / (12.0 * h2);

    h2 = hygrid * hygrid * ct.yside * ct.yside;
    t0 -= 30.0 / (12.0 * h2);
    t1y = 16.0 / (12.0 * h2);
    t2y = -1.0 / (12.0 * h2);

    h2 = hzgrid * hzgrid * ct.zside * ct.zside;
    t0 -= 30.0 / (12.0 * h2);
    t1z = 16.0 / (12.0 * h2);
    t2z = -1.0 / (12.0 * h2);



    for (ix = 2; ix < dimx + 2; ix++)
    {

        for (iy = 2; iy < dimy + 2; iy++)
        {

            for (iz = 2; iz < dimz + 2; iz++)
            {

                work[(ix - 2) * ix1 + (iy - 2) * iy1 + iz - 2] =
                    t0 * dum2[ix * ixs + iy * iys + iz] +
                    t1x * dum2[(ix - 1) * ixs + iy * iys + iz] +
                    t1x * dum2[(ix + 1) * ixs + iy * iys + iz] +
                    t1y * dum2[ix * ixs + (iy - 1) * iys + iz] +
                    t1y * dum2[ix * ixs + (iy + 1) * iys + iz] +
                    t1z * dum2[ix * ixs + iy * iys + iz - 1] +
                    t1z * dum2[ix * ixs + iy * iys + iz + 1] +
                    t2x * dum2[(ix - 2) * ixs + iy * iys + iz] +
                    t2x * dum2[(ix + 2) * ixs + iy * iys + iz] +
                    t2y * dum2[ix * ixs + (iy - 2) * iys + iz] +
                    t2y * dum2[ix * ixs + (iy + 2) * iys + iz] +
                    t2z * dum2[ix * ixs + iy * iys + iz - 2] +
                    t2z * dum2[ix * ixs + iy * iys + iz + 2];

            }                   /* end for */
        }                       /* end for */
    }                           /* end for */

#if MD_TIMERS
    time2 = my_crtc();
    rmg_timings(APPCIL_TIME, (time2 - time1), 0);
#endif

    my_free(dum2);

}                               /* end app6_del2 */


/******/
