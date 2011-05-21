/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/app10_del2f.c *****
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
 *   void app6_del2f(FS0_GRID *f, FP0_GRID *work)
 *   Apply the sixth order Laplacian operator on a orthorhombic grid
 * INPUTS
 *   S0_GRID *f:  real array of (FPX0_GRID+2) * (FPY0_GRID+2) * ((FPZ0_GRID+2)
 *   see main.h for defination of unions S0_GRID and P0_GRID
 * OUTPUT
 *   P0_GRID *work: real array of FPX0_GRID * FPY0_GRID * FPZ0_GRID
 * PARENTS
 *   get_ke.c xcgga.c
 * CHILDREN
 *   trade_image2.c
 * SOURCE
 */



#include "main.h"
#include <float.h>
#include <math.h>




void app10_del2f (FS0_GRID * f, FP0_GRID * work)
{

    int iz, ix, iy;
    REAL h2, t0, t1x, t2x, t3x, t4x;
    REAL t1y, t2y, t3y, t4y;
    REAL t1z, t2z, t3z, t4z;
    REAL *dum2;
    int ixs, iys;
    int ix1, iy1;
    ixs = (FPY0_GRID + 8) * (FPZ0_GRID + 8);
    iys = (FPZ0_GRID + 8);
    ix1 = FPY0_GRID * FPZ0_GRID;
    iy1 = FPZ0_GRID;


    my_malloc (dum2, (FPX0_GRID + 8) * (FPY0_GRID + 8) * (FPZ0_GRID + 8), REAL);

    /*trade_images2f(f, dum2); */
    trade_imagesx (f->s2, dum2, FPX0_GRID, FPY0_GRID, FPZ0_GRID, 4);

    h2 = ct.hxxgrid * ct.hxxgrid * ct.xside * ct.xside;
    t0 = -205.0 / (72.0 * h2);
    t1x = 8.0 / (5.0 * h2);
    t2x = -1.0 / (5.0 * h2);
    t3x = 8.0 / (315.0 * h2);
    t4x = -1.0 / (560.0 * h2);

    h2 = ct.hyygrid * ct.hyygrid * ct.yside * ct.yside;
    t0 -= 205.0 / (72.0 * h2);
    t1y = 8.0 / (5.0 * h2);
    t2y = -1.0 / (5.0 * h2);
    t3y = 8.0 / (315.0 * h2);
    t4y = -1.0 / (560.0 * h2);

    h2 = ct.hzzgrid * ct.hzzgrid * ct.zside * ct.zside;
    t0 -= 205.0 / (72.0 * h2);
    t1z = 8.0 / (5.0 * h2);
    t2z = -1.0 / (5.0 * h2);
    t3z = 8.0 / (315.0 * h2);
    t4z = -1.0 / (560.0 * h2);



    for (ix = 4; ix < FPX0_GRID + 4; ix++)
    {

	    for (iy = 4; iy < FPY0_GRID + 4; iy++)
        {

            for (iz = 4; iz < FPZ0_GRID + 4; iz++)
            {

                work->s1.b[ix - 4][iy - 4][iz - 4] =
                    t0 * dum2[ix * ixs + iy * iys + iz]      +
                    t1x * dum2[(ix - 1) * ixs + iy * iys + iz] +
                    t1x * dum2[(ix + 1) * ixs + iy * iys + iz] +
                    t1y * dum2[ix * ixs + (iy - 1) * iys + iz] +
                    t1y * dum2[ix * ixs + (iy + 1) * iys + iz] +
                    t1z * dum2[ix * ixs + iy * iys + iz - 1]   +
                    t1z * dum2[ix * ixs + iy * iys + iz + 1]   +
                    t2x * dum2[(ix - 2) * ixs + iy * iys + iz] +
                    t2x * dum2[(ix + 2) * ixs + iy * iys + iz] +
                    t2y * dum2[ix * ixs + (iy - 2) * iys + iz] +
                    t2y * dum2[ix * ixs + (iy + 2) * iys + iz] +
                    t2z * dum2[ix * ixs + iy * iys + iz - 2]   +
                    t2z * dum2[ix * ixs + iy * iys + iz + 2]   +
                    t3x * dum2[(ix - 3) * ixs + iy * iys + iz] +
                    t3x * dum2[(ix + 3) * ixs + iy * iys + iz] +
                    t3y * dum2[ix * ixs + (iy - 3) * iys + iz] +
                    t3y * dum2[ix * ixs + (iy + 3) * iys + iz] +
                    t3z * dum2[ix * ixs + iy * iys + iz - 3]   +
                    t3z * dum2[ix * ixs + iy * iys + iz + 3]   +
                    t4x * dum2[(ix - 4) * ixs + iy * iys + iz] +
                    t4x * dum2[(ix + 4) * ixs + iy * iys + iz] +
                    t4y * dum2[ix * ixs + (iy - 4) * iys + iz] +
                    t4y * dum2[ix * ixs + (iy + 4) * iys + iz] +
                    t4z * dum2[ix * ixs + iy * iys + iz - 4]   +
                    t4z * dum2[ix * ixs + iy * iys + iz + 4];   


            }                   /* end for */
        }                       /* end for */
    }                           /* end for */

    my_free (dum2);



}                               /* end app6_del2 */


/******/
