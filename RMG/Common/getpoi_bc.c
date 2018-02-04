/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/




#include "main.h"
#include "common_prototypes.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



void getpoi_bc (double * rho, double * vh_bc, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, idx, stop, PX0_GRID, PY0_GRID, PZ0_GRID;
    int ix1, iy1;
    int ixdim, iydim, izdim;
    int pex, pey, pez;
    double ir2, q, px, py, pz, sxx, syy, szz, sxy, syz, szx, temp;
    double xc, yc, zc, x, y, z;
    double ax[3], bx[3];
    double xoff, yoff, zoff;
    double *mask;
    double hxgrid, hygrid, hzgrid;
    int incx, incy, incz;

    PX0_GRID = get_PX0_GRID();
    PY0_GRID = get_PY0_GRID();
    PZ0_GRID = get_PZ0_GRID();

    hxgrid = get_hxgrid();
    hygrid = get_hygrid();
    hzgrid = get_hzgrid();

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
    my_malloc (mask, stop, double);

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


    pe2xyz (pct.gridpe, &pex, &pey, &pez);

    xc = pex * hxgrid * PX0_GRID;

    for (ix = 0; ix < PX0_GRID; ix++)
    {

        yc = pey * hygrid * PY0_GRID;
        for (iy = 0; iy < PY0_GRID; iy++)
        {

            zc = pez * hzgrid * PZ0_GRID;
            for (iz = 0; iz < PZ0_GRID; iz++)
            {

                ax[0] = xc - 0.5;
                ax[1] = yc - 0.5;
                ax[2] = zc - 0.5;
                to_cartesian (ax, bx);
                x = bx[0];
                y = bx[1];
                z = bx[2];

                temp = rho[ix * incx + iy * incy + iz] * get_vel();
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

                zc += hzgrid;

            }                   /* end for */

            yc += hygrid;

        }                       /* end for */

        xc += hxgrid;

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
    q = real_sum_all (q, pct.grid_comm);
    px = real_sum_all (px, pct.grid_comm);
    py = real_sum_all (py, pct.grid_comm);
    pz = real_sum_all (pz, pct.grid_comm);

    sxx = real_sum_all (sxx, pct.grid_comm);
    syy = real_sum_all (syy, pct.grid_comm);
    szz = real_sum_all (szz, pct.grid_comm);

    sxy = real_sum_all (sxy, pct.grid_comm);
    syz = real_sum_all (syz, pct.grid_comm);
    szx = real_sum_all (szx, pct.grid_comm);
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
    xc = pex * hxgrid * ixdim - xoff - hxgrid;
    for (ix = 0; ix < ixdim + 2; ix++)
    {

        yc = pey * hygrid * iydim - yoff - hygrid;
        for (iy = 0; iy < iydim + 2; iy++)
        {

            for (iz = 0; iz < izdim + 2; iz += incz)
            {

                zc = pez * hzgrid * izdim + ((double) iz) * hzgrid - zoff - hzgrid;
                ax[0] = xc;
                ax[1] = yc;
                ax[2] = zc;
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

                            ax[0] = xc + (double) ix1;
                            ax[1] = yc + (double) iy1;
                            ax[2] = zc;
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

            yc += hygrid;

        }                       /* end for */

        xc += hxgrid;

    }                           /* end for */

    for (idx = 0; idx < stop; idx++)
    {
        vh_bc[idx] = vh_bc[idx] * mask[idx];
    }

    my_free (mask);

}                               /* end getpoi_bc.c */

/******/
