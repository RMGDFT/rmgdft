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

    double sumx, sumy, sumz;

    int ix, iy, iz, ixx, iyy, izz;
    int xstart, ystart, zstart, xend, yend, zend;
    int ion, idx;
    int ilow, jlow, klow, ihi, jhi, khi;
    int dimx, dimy, dimz;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID, FP0_BASIS;
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
    FP0_BASIS = FPX0_GRID * FPY0_GRID * FPZ0_GRID;
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


    double *force_nlcc, *gx, *gy, *gz;
    double rho_at_r;
    force_nlcc = (double *)malloc(ct.num_ions *3 * sizeof(double));
    gx = (double *)malloc(FP0_BASIS * sizeof(double));
    gy = (double *)malloc(FP0_BASIS * sizeof(double));
    gz = (double *)malloc(FP0_BASIS * sizeof(double));


    app_grad (vxc, gx, gy, gz, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid);

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
                                    rho_at_r = AtomicInterpolate (&sp->rhocorelig[0], r);
                                        sumx +=  rho_at_r * gx[idx]; 
                                    sumy +=  rho_at_r * gy[idx]; 
                                    sumz +=  rho_at_r * gz[idx]; 


                                }   /* end if */
                            }       /* end for */
                        }           /* end for */
                    }               /* end for */
                }                   /* end if */
            }


            force_nlcc[ion * 3 + 0] = -sumx * get_vel_f();
            force_nlcc[ion * 3 + 1] = -sumy * get_vel_f();
            force_nlcc[ion * 3 + 2] = -sumz * get_vel_f();


        }                       /* end if */


    }                           /* end for */

    idx = 3* ion;
    global_sums (force_nlcc, &idx, pct.grid_comm);
    double fac_spin = 1.0/(1.0 + ct.spin_flag);
    /* factor 0.5 is because when calculating exchange correlation
       half of nonlinear core corection charge is added to spin up and down density */

    for (ion = 0; ion < ct.num_ions; ion++)
    {
    //    printf("\n force_nlcc  %d %e %e %e", ion, force_nlcc[ion*3], force_nlcc[ion*3+1], force_nlcc[ion*3+2]);

        /* Generate ion pointer */
        iptr = &ct.ions[ion];
        iptr->force[ct.fpt[0]][0] += force_nlcc[ion*3 + 0] * fac_spin;
        iptr->force[ct.fpt[0]][1] += force_nlcc[ion*3 + 1] * fac_spin;
        iptr->force[ct.fpt[0]][2] += force_nlcc[ion*3 + 2] * fac_spin;
    }

    free(force_nlcc);
    free(gx);
    free(gy);
    free(gz);

}                               /* end nlccforce */

/******/




