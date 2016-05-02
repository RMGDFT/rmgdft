/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"



void lforce (double * rho, double * vh, double *force)
{

    int ix, iy, iz, ixx, iyy, izz;
    int xstart, ystart, zstart, xend, yend, zend;
    int ion, idx;
    int ilow, jlow, klow, ihi, jhi, khi;
    int dimx, dimy, dimz;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID;
    int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
    int FNX_GRID, FNY_GRID, FNZ_GRID;

    double r, Zv, rc, rc2, rcnorm, t1;
    double x[3];
    double hxxgrid, hyygrid, hzzgrid;
    double xside, yside, zside;
    SPECIES *sp;
    ION *iptr;

    double bx[3], norm1, fx, fy, fz;


    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();
    xside = get_xside();
    yside = get_yside();
    zside = get_zside();

    FPX0_GRID = get_FPX0_GRID();
    FPY0_GRID = get_FPY0_GRID();
    FPZ0_GRID = get_FPZ0_GRID();
    FPX_OFFSET = get_FPX_OFFSET();
    FPY_OFFSET = get_FPY_OFFSET();
    FPZ_OFFSET = get_FPZ_OFFSET();
    FNX_GRID = get_FNX_GRID();
    FNY_GRID = get_FNY_GRID();
    FNZ_GRID = get_FNZ_GRID();


    ilow = FPX_OFFSET;
    jlow = FPY_OFFSET;
    klow = FPZ_OFFSET;
    ihi = ilow + FPX0_GRID;
    jhi = jlow + FPY0_GRID;
    khi = klow + FPZ0_GRID;

    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Generate ion pointer */
        iptr = &ct.ions[ion];


        /* Get species type */
        sp = &ct.sp[iptr->species];


        Zv = sp->zvalence;
        rc = sp->rc;
        rc2 = rc * rc;
        rcnorm = rc * rc * rc * pow (PI, 1.5);
        rcnorm = ONE / rcnorm;
        rc2 = sp->rc * sp->rc;
        norm1 = -Zv * rcnorm / rc2;


        dimx =  sp->lradius/(hxxgrid*xside);
        dimy =  sp->lradius/(hyygrid*yside);
        dimz =  sp->lradius/(hzzgrid*zside);

        dimx = dimx * 2 + 1;
        dimy = dimy * 2 + 1;
        dimz = dimz * 2 + 1;


        xstart = iptr->xtal[0] / hxxgrid - dimx/2;
        xend = xstart + dimx;
        ystart = iptr->xtal[1] / hyygrid - dimy/2;
        yend = ystart + dimy;
        zstart = iptr->xtal[2] / hzzgrid - dimz/2;
        zend = zstart + dimz;

        fx = fy = fz = 0.0;

        for (ix = xstart; ix < xend; ix++)
        {
            // fold the grid into the unit cell
            ixx = (ix + 20 * FNX_GRID) % FNX_GRID;
            if(ixx >= ilow && ixx < ihi)
            {

                for (iy = ystart; iy < yend; iy++)
                {
                    // fold the grid into the unit cell
                    iyy = (iy + 20 * FNY_GRID) % FNY_GRID;
                    if(iyy >= jlow && iyy < jhi)
                    {
                        for (iz = zstart; iz < zend; iz++)
                        {
                            // fold the grid into the unit cell
                            izz = (iz + 20 * FNZ_GRID) % FNZ_GRID;
                            if(izz >= klow && izz < khi)
                            {

                                idx = (ixx-ilow) * FPY0_GRID * FPZ0_GRID + (iyy-jlow) * FPZ0_GRID + izz-klow;
                                x[0] = ix * hxxgrid - iptr->xtal[0];
                                x[1] = iy * hyygrid - iptr->xtal[1];
                                x[2] = iz * hzzgrid - iptr->xtal[2];
                                r = metric (x);





                                to_cartesian (x, bx);
                                r = metric (x);

                                t1 = 2.0 * norm1 * exp (-r * r / rc2);
                                fx += bx[0] * t1 * vh[idx];
                                fy += bx[1] * t1 * vh[idx];
                                fz += bx[2] * t1 * vh[idx];

                                t1 = AtomicInterpolate (&sp->drlocalig[0], r);
                                fx += -t1 * bx[0] / (r+1.0e-10) * rho[idx];
                                fy += -t1 * bx[1] / (r+1.0e-10) * rho[idx];
                                fz += -t1 * bx[2] / (r+1.0e-10) * rho[idx];


                            }
                        }
                    }
                }
            }

        }                       /* end for */

        force[ion * 3 + 0] = -get_vel_f() * fx;
        force[ion * 3 + 1] = -get_vel_f() * fy;
        force[ion * 3 + 2] = -get_vel_f() * fz;

    }                           /* end for */




}                               /* end lforce */

/******/


