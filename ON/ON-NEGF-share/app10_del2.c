/************************** SVN Revision Information **************************
 **    $Id  $    **
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


void app10_del2(rmg_double_t * f, rmg_double_t * work, int dimx, int dimy, int dimz,
                rmg_double_t hxgrid, rmg_double_t hygrid, rmg_double_t hzgrid)
{

    int iz, ix, iy;
    rmg_double_t h2, t0, t1x, t2x, t3x, t4x;
    rmg_double_t t1y, t2y, t3y, t4y;
    rmg_double_t t1z, t2z, t3z, t4z;
    rmg_double_t *dum2;
    int ixs, iys, ix1, iy1;

#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc();
#endif

    ixs = (dimy + 8) * (dimz + 8);
    iys = (dimz + 8);
    ix1 = dimy * dimz;
    iy1 = dimz;

    my_malloc_init( dum2, (dimx + 8) * (dimy + 8) * (dimz + 8), rmg_double_t );
    fill_orbit_borders4(dum2, f, dimx, dimy, dimz);

    h2 = hxgrid * hxgrid * get_xside() * get_xside();
    t0 = -205.0 / (72.0 * h2);
    t1x = 8.0 / (5.0 * h2);
    t2x = -1.0 / (5.0 * h2);
    t3x = 8.0 / (315.0 * h2);
    t4x = -1.0 / (560.0 * h2);

    h2 = hygrid * hygrid * get_yside() * get_yside();
    t0 -= 205.0 / (72.0 * h2);
    t1y = 8.0 / (5.0 * h2);
    t2y = -1.0 / (5.0 * h2);
    t3y = 8.0 / (315.0 * h2);
    t4y = -1.0 / (560.0 * h2);

    h2 = hzgrid * hzgrid * get_zside() * get_zside();
    t0 -= 205.0 / (72.0 * h2);
    t1z = 8.0 / (5.0 * h2);
    t2z = -1.0 / (5.0 * h2);
    t3z = 8.0 / (315.0 * h2);
    t4z = -1.0 / (560.0 * h2);



    for (ix = 4; ix < dimx + 4; ix++)
    {

        for (iy = 4; iy < dimy + 4; iy++)
        {

            for (iz = 4; iz < dimz + 4; iz++)
            {

                work[(ix - 4) * ix1 + (iy - 4) * iy1 + iz - 4] =
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
                    t2z * dum2[ix * ixs + iy * iys + iz + 2] +
                    t3x * dum2[(ix - 3) * ixs + iy * iys + iz] +
                    t3x * dum2[(ix + 3) * ixs + iy * iys + iz] +
                    t3y * dum2[ix * ixs + (iy - 3) * iys + iz] +
                    t3y * dum2[ix * ixs + (iy + 3) * iys + iz] +
                    t3z * dum2[ix * ixs + iy * iys + iz - 3] +
                    t3z * dum2[ix * ixs + iy * iys + iz + 3] +
                    t4x * dum2[(ix - 4) * ixs + iy * iys + iz] +
                    t4x * dum2[(ix + 4) * ixs + iy * iys + iz] +
                    t4y * dum2[ix * ixs + (iy - 4) * iys + iz] +
                    t4y * dum2[ix * ixs + (iy + 4) * iys + iz] +
                    t4z * dum2[ix * ixs + iy * iys + iz - 4] +
                    t4z * dum2[ix * ixs + iy * iys + iz + 4];

            }                   /* end for */
        }                       /* end for */
    }                           /* end for */

#if MD_TIMERS
    time2 = my_crtc();
    rmg_timings(APPCIL_TIME, (time2 - time1));
#endif

    my_free(dum2);

}                               /* end app8_del2 */


/******/
