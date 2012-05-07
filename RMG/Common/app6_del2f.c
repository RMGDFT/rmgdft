/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/app6_del2f.c *****
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




void app6_del2f (REAL *rho, FP0_GRID * work)
{

    int iz, ix, iy;
    REAL h2, t0, t1x, t2x;
    REAL t1y, t2y;
    REAL t1z, t2z;
    FSS0_GRID *f;


    my_malloc (f, 1, FSS0_GRID);

    trade_imagesx (rho, &f->b[0][0][0], FPX0_GRID, FPY0_GRID, FPZ0_GRID, 2, CENTRAL_FD);

    h2 = ct.hxxgrid * ct.hxxgrid * ct.xside * ct.xside;
    t0 = -30.0 / (12.0 * h2);
    t1x = 16.0 / (12.0 * h2);
    t2x = -1.0 / (12.0 * h2);

    h2 = ct.hyygrid * ct.hyygrid * ct.yside * ct.yside;
    t0 -= 30.0 / (12.0 * h2);
    t1y = 16.0 / (12.0 * h2);
    t2y = -1.0 / (12.0 * h2);

    h2 = ct.hzzgrid * ct.hzzgrid * ct.zside * ct.zside;
    t0 -= 30.0 / (12.0 * h2);
    t1z = 16.0 / (12.0 * h2);
    t2z = -1.0 / (12.0 * h2);



    for (ix = 2; ix < FPX0_GRID + 2; ix++)
    {

        for (iy = 2; iy < FPY0_GRID + 2; iy++)
        {

            for (iz = 2; iz < FPZ0_GRID + 2; iz++)
            {

                work->s1.b[ix - 2][iy - 2][iz - 2] = t0 * f->b[ix][iy][iz] +
                    t1x * f->b[ix - 1][iy][iz] +
                    t1x * f->b[ix + 1][iy][iz] +
                    t1y * f->b[ix][iy - 1][iz] +
                    t1y * f->b[ix][iy + 1][iz] +
                    t1z * f->b[ix][iy][iz - 1] +
                    t1z * f->b[ix][iy][iz + 1] +
                    t2x * f->b[ix - 2][iy][iz] +
                    t2x * f->b[ix + 2][iy][iz] +
                    t2y * f->b[ix][iy - 2][iz] +
                    t2y * f->b[ix][iy + 2][iz] +
                    t2z * f->b[ix][iy][iz - 2] + t2z * f->b[ix][iy][iz + 2];

            }                   /* end for */
        }                       /* end for */
    }                           /* end for */

    my_free(f);

}                               /* end app6_del2 */


/******/
