/************************** SVN Revision Information **************************
 **    $Id: mg_restrict.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
/****f* QMD-MGDFT/mg_restrict.c *****
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
 *   void mg_restrict(REAL *full, REAL *half, int dimx, int dimy, int dimz)
 *   get half dimensioned array in the corse grid from full dimensioned 
 *   array in the fine grid. The returned values are smoothed.
 * INPUTS
 *   full[(dimx+2)*(dimy+2)*(dimz+2)]: array in the fine grid
 *      image value must be there
 *   dimx, dimy, dimz: dimensions of array in the fine grid
 * OUTPUT
 *   half[(dimx/2+2)*(dimy/2+2)*(dimz/2+2)] array in the corse grid
 * PARENTS
 *   mgrid_solv.c
 * CHILDREN
 *   nothing 
 * SOURCE
 */



#include "md.h"
#include <float.h>
#include <math.h>





void mg_restrict(REAL * full, REAL * half, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int incz, incy, incx, incz2, incy2, incx2;
    int x0, xp, xm, y0, yp, ym, z0, zp, zm;
    REAL scale, face, corner, edge;



    incz = 1;
    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);

    incz2 = 1;
    incy2 = dimz / 2 + 2;
    incx2 = (dimz / 2 + 2) * (dimy / 2 + 2);


    switch (ct.ibrav)
    {

    case CUBIC_PRIMITIVE:
    case CUBIC_FC:
    case ORTHORHOMBIC_PRIMITIVE:
    case HEXAGONAL:

        scale = ONE / 64.0;
        for (ix = 1; ix <= dimx / 2; ix++)
        {

            x0 = 2 * ix - 1;
            xp = x0 + 1;
            xm = x0 - 1;

            for (iy = 1; iy <= dimy / 2; iy++)
            {

                y0 = 2 * iy - 1;
                yp = y0 + 1;
                ym = y0 - 1;

                for (iz = 1; iz <= dimz / 2; iz++)
                {

                    z0 = 2 * iz - 1;
                    zp = z0 + 1;
                    zm = z0 - 1;

                    face = full[xm * incx + y0 * incy + z0] +
                        full[xp * incx + y0 * incy + z0] +
                        full[x0 * incx + ym * incy + z0] +
                        full[x0 * incx + yp * incy + z0] +
                        full[x0 * incx + y0 * incy + zm] + full[x0 * incx + y0 * incy + zp];

                    corner =
                        full[xm * incx + ym * incy + zm] +
                        full[xm * incx + ym * incy + zp] +
                        full[xm * incx + yp * incy + zm] +
                        full[xm * incx + yp * incy + zp] +
                        full[xp * incx + ym * incy + zm] +
                        full[xp * incx + ym * incy + zp] +
                        full[xp * incx + yp * incy + zm] + full[xp * incx + yp * incy + zp];

                    edge = full[xm * incx + y0 * incy + zm] +
                        full[xm * incx + ym * incy + z0] +
                        full[xm * incx + yp * incy + z0] +
                        full[xm * incx + y0 * incy + zp] +
                        full[x0 * incx + ym * incy + zm] +
                        full[x0 * incx + yp * incy + zm] +
                        full[x0 * incx + ym * incy + zp] +
                        full[x0 * incx + yp * incy + zp] +
                        full[xp * incx + y0 * incy + zm] +
                        full[xp * incx + ym * incy + z0] +
                        full[xp * incx + yp * incy + z0] + full[xp * incx + y0 * incy + zp];


                    half[ix * incx2 + iy * incy2 + iz] =
                        scale * (8.0 * full[x0 * incx + y0 * incy + z0] +
                                 4.0 * face + 2.0 * edge + corner);


                }               /* end for */

            }                   /* end for */

        }                       /* end for */

        break;


    case CUBIC_BC:

        scale = ONE / 52.0;

        for (ix = 1; ix <= dimx / 2; ix++)
        {

            x0 = 2 * ix - 1;
            xp = x0 + 1;
            xm = x0 - 1;

            for (iy = 1; iy <= dimy / 2; iy++)
            {

                y0 = 2 * iy - 1;
                yp = y0 + 1;
                ym = y0 - 1;

                for (iz = 1; iz <= dimz / 2; iz++)
                {

                    z0 = 2 * iz - 1;
                    zp = z0 + 1;
                    zm = z0 - 1;

                    face = full[xm * incx + ym * incy + z0] +
                        full[xm * incx + y0 * incy + zm] +
                        full[x0 * incx + ym * incy + zm] +
                        full[x0 * incx + yp * incy + zp] +
                        full[xp * incx + y0 * incy + zp] + full[xp * incx + yp * incy + z0];

                    corner =
                        full[xm * incx + ym * incy + zm] +
                        full[xm * incx + y0 * incy + z0] +
                        full[x0 * incx + ym * incy + z0] +
                        full[x0 * incx + y0 * incy + zm] +
                        full[x0 * incx + y0 * incy + zp] +
                        full[x0 * incx + yp * incy + z0] +
                        full[xp * incx + y0 * incy + z0] + full[xp * incx + yp * incy + zp];


                    half[ix * incx2 + iy * incy2 + iz] =
                        scale * (8.0 * full[x0 * incx + y0 * incy + z0] +
                                 4.0 * corner + 2.0 * face);


                }               /* end for */

            }                   /* end for */

        }                       /* end for */

        break;

    case 20:

        scale = ONE / 80.0;
        for (ix = 1; ix <= dimx / 2; ix++)
        {

            x0 = 2 * ix - 1;
            xp = x0 + 1;
            xm = x0 - 1;

            for (iy = 1; iy <= dimy / 2; iy++)
            {

                y0 = 2 * iy - 1;
                yp = y0 + 1;
                ym = y0 - 1;

                for (iz = 1; iz <= dimz / 2; iz++)
                {

                    z0 = 2 * iz - 1;
                    zp = z0 + 1;
                    zm = z0 - 1;

                    face = full[xm * incx + ym * incy + z0] +
                        full[xm * incx + y0 * incy + zm] +
                        full[x0 * incx + ym * incy + zm] +
                        full[x0 * incx + yp * incy + zp] +
                        full[xp * incx + y0 * incy + zp] + full[xp * incx + yp * incy + z0];

                    edge =
                        full[xm * incx + y0 * incy + z0] +
                        full[xm * incx + y0 * incy + zp] +
                        full[xm * incx + yp * incy + z0] +
                        full[x0 * incx + ym * incy + z0] +
                        full[x0 * incx + ym * incy + zp] +
                        full[x0 * incx + y0 * incy + zm] +
                        full[x0 * incx + y0 * incy + zp] +
                        full[x0 * incx + yp * incy + zm] +
                        full[x0 * incx + yp * incy + z0] +
                        full[xp * incx + ym * incy + z0] +
                        full[xp * incx + y0 * incy + zm] + full[xp * incx + y0 * incy + z0];


                    half[ix * incx2 + iy * incy2 + iz] =
                        scale * (8.0 * full[x0 * incx + y0 * incy + z0] + 5.0 * edge + 2.0 * face);


                }               /* end for */

            }                   /* end for */

        }                       /* end for */

        break;

    default:
        error_handler("Lattice type not programmed");

    }                           /* end switch */


}                               /* end mg_restrict */


/******/
