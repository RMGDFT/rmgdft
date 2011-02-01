/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/getpoi_bc.c *****
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
 *   void getpoi_bc(REAL *rho, REAL *vh_bc, int dimx, int dimy, int dimz)
 *   Computes the hartree potential boundary conditions via multipole
 *   expansion about the center when using cluster or surface boundary
 *   conditions
 * INPUTS
 *   rho: total charge density
 *   dimx, dimy, dimz: dimensions of rho 
 * OUTPUT
 *   vh_bc:  Hartree potential
 * PARENTS
 *   get_vh.c
 * CHILDREN
 *   set_bc.c pe2xyz.c
 * SOURCE
 */




#include "main.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



void getpoi_bc (REAL * rho, REAL * vh_bc, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, idx, stop;
    int ix1, iy1;
    int ixdim, iydim, izdim;
    int pex, pey, pez;
    REAL ir2, q, px, py, pz, sxx, syy, szz, sxy, syz, szx, temp;
    REAL r, xc, yc, zc, x, y, z;
    REAL ax[3], bx[3];
    REAL xoff, yoff, zoff;
    REAL *mask;
    int incx, incy, incz;

    ixdim = 2 * dimx;
    iydim = 2 * dimy;
    izdim = 2 * dimz;
    xoff = 0.5;
    yoff = 0.5;
    zoff = 0.5;
    if (ct.boundaryflag == SURFACE)
    {

        ixdim = dimx;
        iydim = dimy;
        xoff = 0.0;
        yoff = 0.0;

    }                           /* end if */

    stop = (ixdim + 2) * (iydim + 2) * (izdim + 2);
    my_malloc (mask, stop, REAL);

    for (idx = 0; idx < stop; idx++)
        mask[idx] = 0.0;

    set_bc (mask, ixdim, iydim, izdim, 1, 1.0);


    incy = dimz;
    incx = dimy * dimz;


    /*  compute integral of various moments of r and rho  */
    q = ZERO;
    px = ZERO;
    py = ZERO;
    pz = ZERO;
    sxx = ZERO;
    syy = ZERO;
    szz = ZERO;
    sxy = ZERO;
    syz = ZERO;
    szx = ZERO;


    pe2xyz (pct.thispe, &pex, &pey, &pez);

    xc = pex * ct.hxgrid * PX0_GRID;

    for (ix = 0; ix < PX0_GRID; ix++)
    {

        yc = pey * ct.hygrid * PY0_GRID;
        for (iy = 0; iy < PY0_GRID; iy++)
        {

            zc = pez * ct.hzgrid * PZ0_GRID;
            for (iz = 0; iz < PZ0_GRID; iz++)
            {

                ax[0] = xc - 0.5;
                ax[1] = yc - 0.5;
                ax[2] = zc - 0.5;
                r = metric (ax);
                to_cartesian (ax, bx);
                x = bx[0];
                y = bx[1];
                z = bx[2];

                temp = rho[ix * incx + iy * incy + iz] * ct.vel;
                q += temp;

                px = px + x * temp;
                py = py + y * temp;
                pz = pz + z * temp;
                sxx = sxx + (2.0 * x * x - (y * y + z * z)) * temp;
                syy = syy + (2.0 * y * y - (x * x + z * z)) * temp;
                szz = szz + (2.0 * z * z - (y * y + x * x)) * temp;
                sxy = sxy + x * y * temp;
                syz = syz + y * z * temp;
                szx = szx + z * x * temp;

                zc += ct.hzgrid;

            }                   /* end for */

            yc += ct.hygrid;

        }                       /* end for */

        xc += ct.hxgrid;

    }                           /* end for */



    /* some correction factors for the off-diagonal quadrapole moments  */
    /* (Note: that these differ from Jackson's definitions by 1/2)      */
    sxx = sxx / 2.0;
    syy = syy / 2.0;
    szz = szz / 2.0;
    sxy *= 1.5;
    syz *= 1.5;
    szx *= 1.5;


    /* Sum these up over all processors */
    q = real_sum_all (q);
    px = real_sum_all (px);
    py = real_sum_all (py);
    pz = real_sum_all (pz);

    sxx = real_sum_all (sxx);
    syy = real_sum_all (syy);
    szz = real_sum_all (szz);

    sxy = real_sum_all (sxy);
    syz = real_sum_all (syz);
    szx = real_sum_all (szx);
    printf ("QQ = %12.6f\n", q);
    printf ("PX = %12.6f\n", px);
    printf ("PY = %12.6f\n", py);
    printf ("PZ = %12.6f\n", pz);


    incy = izdim + 2;
    incx = (izdim + 2) * (iydim + 2);
    incz = 1;
    if (ct.boundaryflag == SURFACE)
        incz = izdim + 1;

    /* Set the boundary condition on the surface of the grid */
    xc = pex * ct.hxgrid * ixdim - xoff - ct.hxgrid;
    for (ix = 0; ix < ixdim + 2; ix++)
    {

        yc = pey * ct.hygrid * iydim - yoff - ct.hygrid;
        for (iy = 0; iy < iydim + 2; iy++)
        {

            for (iz = 0; iz < izdim + 2; iz += incz)
            {

                zc = pez * ct.hzgrid * izdim + ((REAL) iz) * ct.hzgrid - zoff - ct.hzgrid;
                ax[0] = xc;
                ax[1] = yc;
                ax[2] = zc;
                r = metric (ax);
                to_cartesian (ax, bx);
                x = bx[0];
                y = bx[1];
                z = bx[2];

                ir2 = 1.0 / (x * x + y * y + z * z + 1.0e-10);

                vh_bc[ix * incx + iy * incy + iz] = sqrt (ir2) *
                    (q +
                     ir2 * (x * px + y * py + z * pz +
                            ir2 * (x * (sxx * x + sxy * y + szx * z) +
                                   y * (sxy * x + syy * y + syz * z) + z * (szx * x + syz * y +
                                                                            szz * z))));


                /* If we have surface boundary conditions we need to add 
                 * the contributions from adjoining cells in the x and y planes */
                if (ct.boundaryflag == SURFACE)
                {

                    vh_bc[ix * incx + iy * incy + iz] = 0.0;
                    for (ix1 = -2; ix1 <= 2; ix1++)
                    {
                        for (iy1 = -2; iy1 <= 2; iy1++)
                        {

                            ax[0] = xc + (REAL) ix1;
                            ax[1] = yc + (REAL) iy1;
                            ax[2] = zc;
                            r = metric (ax);
                            to_cartesian (ax, bx);
                            x = bx[0];
                            y = bx[1];
                            z = bx[2];
                            ir2 = 1.0 / (x * x + y * y + z * z + 1.0e-10);
                            vh_bc[ix * incx + iy * incy + iz] += sqrt (ir2) *
                                (q +
                                 ir2 * (x * px + y * py + z * pz +
                                        ir2 * (x * (sxx * x + sxy * y + szx * z) +
                                               y * (sxy * x + syy * y + syz * z) + z * (szx * x +
                                                                                        syz * y +
                                                                                        szz * z))));

                        }       /* end for */

                    }           /* end for */

                }               /* end if */

            }                   /* end for */

            yc += ct.hygrid;

        }                       /* end for */

        xc += ct.hxgrid;

    }                           /* end for */

    for (idx = 0; idx < stop; idx++)
    {
        vh_bc[idx] = vh_bc[idx] * mask[idx];
    }

    my_free (mask);

}                               /* end getpoi_bc.c */

/******/
