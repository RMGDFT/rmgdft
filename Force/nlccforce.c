/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/nlccforce.c *****
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
 *   void nlccforce(double *rho, double *vxc)
 *   Evaluates the ionic force component due to the non-linear core
 *   correction terms. 
 * INPUTS
 *   rho: total charge density
 *   vxc: exchange correlation potential 
 *   core charge density is obtained from ct.sp->rhocorelig
 * OUTPUT
 *   forces are added to structure ct.ion
 * PARENTS
 *   force.c
 * CHILDREN
 *   get_index.c to_crystal.c to_cartesian.c linint.c
 * SOURCE
 */




#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"



void nlccforce (double * rho, double * vxc)
{

    int ishift;
    double axs[3], bx[3];
    double shift[4];
    double fl[3];
    double deltac;
    double sumxc2, sumx, sumy, sumz;

    int ix, iy, iz, ixx, iyy, izz;
    int xstart, ystart, zstart, xend, yend, zend;
    int ion, idx;
    int ilow, jlow, klow, ihi, jhi, khi;
    int dimx, dimy, dimz;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID;
    int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
    int FNX_GRID, FNY_GRID, FNZ_GRID;

    double r;
    double x[3];
    double hxxgrid, hyygrid, hzzgrid;
    double xside, yside, zside;
    SPECIES *sp;
    ION *iptr;

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


    deltac = ct.hmaxgrid / 200.0 / (double) get_FG_RATIO();
    shift[0] = -TWO * deltac;
    shift[1] = TWO * deltac;
    shift[2] = -deltac;
    shift[3] = deltac;

    sumxc2 = TWO * (shift[0] * shift[0] + shift[2] * shift[2]);

    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];

        if (sp->nlccflag)
        {

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

            sumx = sumy = sumz = 0.0;
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


                                    for (ishift = 0; ishift < 4; ishift++)
                                    {
                                        axs[0] = bx[0] - shift[ishift];
                                        axs[1] = bx[1];
                                        axs[2] = bx[2];
                                        r = sqrt (axs[0] * axs[0] + axs[1] * axs[1] + axs[2] * axs[2]);
                                        sumx +=   AtomicInterpolate (&sp->rhocorelig[0], r) * shift[ishift] * vxc[idx];
                                    }       /* end for */


                                    for (ishift = 0; ishift < 4; ishift++)
                                    {
                                        axs[0] = bx[0];
                                        axs[1] = bx[1] - shift[ishift];
                                        axs[2] = bx[2];
                                        r = sqrt (axs[0] * axs[0] + axs[1] * axs[1] + axs[2] * axs[2]);
                                        sumy +=   AtomicInterpolate (&sp->rhocorelig[0], r) * shift[ishift] * vxc[idx];
                                    }       /* end for */

                                    for (ishift = 0; ishift < 4; ishift++)
                                    {
                                        axs[0] = bx[0];
                                        axs[1] = bx[1];
                                        axs[2] = bx[2] - shift[ishift];
                                        r = sqrt (axs[0] * axs[0] + axs[1] * axs[1] + axs[2] * axs[2]);
                                        sumz +=   AtomicInterpolate (&sp->rhocorelig[0], r) * shift[ishift] * vxc[idx];
                                    }       /* end for */
                                }   /* end if */
                            }       /* end for */
                        }           /* end for */
                    }               /* end for */
                }                   /* end if */
            }


            fl[0] = -sumx/sumxc2 * get_vel_f();
            fl[1] = -sumy/sumxc2 * get_vel_f();
            fl[2] = -sumz/sumxc2 * get_vel_f();



            idx = 3;
            global_sums (fl, &idx, pct.grid_comm);
            if (ct.spin_flag)
            {
                fl[0] *= 0.5; 
                fl[1] *= 0.5; 
                fl[2] *= 0.5; 
            }
            /* factor 0.5 is because when calculating exchange correlation
               half of nonlinear core corection charge is added to spin up and down density */


            iptr->force[ct.fpt[0]][0] += fl[0];
            iptr->force[ct.fpt[0]][1] += fl[1];
            iptr->force[ct.fpt[0]][2] += fl[2];


        }                       /* end if */


    }                           /* end for */


}                               /* end nlccforce */

/******/




