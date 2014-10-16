/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"


void init_nuc (rmg_double_t * vnuc_f, rmg_double_t * rhoc_f, rmg_double_t * rhocore_f)
{

    int ix, iy, iz, ixx, iyy, izz;
    int xstart, ystart, zstart, xend, yend, zend;
    int ion, idx;
    int ilow, jlow, klow, ihi, jhi, khi;
    int dimx, dimy, dimz;
    int FP0_BASIS;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID;
    int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
    int FNX_GRID, FNY_GRID, FNZ_GRID;

    rmg_double_t r, Zv, rc, rc2, rcnorm, t1;
    rmg_double_t x[3], invdr;
    rmg_double_t hxxgrid, hyygrid, hzzgrid;
    double xside, yside, zside;
    SPECIES *sp;
    ION *iptr;

    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();
    xside = get_xside();
    yside = get_yside();
    zside = get_zside();

    FP0_BASIS = get_FP0_BASIS();
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


    /* Initialize the compensating charge array and the core charge array */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        rhoc_f[idx] = ct.background_charge / FP0_BASIS / get_vel_f() / NPES;
        rhocore_f[idx] = 0.0;
        vnuc_f[idx] = 0.0;
    }

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
        rcnorm = 1.0 / rcnorm;

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

        invdr = 1.0 / sp->drlig;

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

                                vnuc_f[idx] += linint (&sp->localig[0], r, invdr);
                                rhoc_f[idx] += Zv * exp (-r * r / rc2) * rcnorm;

                                if (sp->nlccflag)
                                    rhocore_f[idx] += linint (&sp->rhocorelig[0], r, invdr);


                            }                           /* end for */

                        }
                    }
                }
            }

        }
    }
    /* Check compensating charges */
    ct.crho = 0.0;
    for (idx = 0; idx < FP0_BASIS; idx++)
        ct.crho += rhoc_f[idx];


    ct.crho = ct.crho * get_vel_f();
    ct.crho = real_sum_all (ct.crho, pct.grid_comm);  /* sum over pct.grid_comm  */

    if (pct.imgpe==0)
        printf("\nCompensating charge is %.8e\n", ct.crho);


    t1 = 0.0;
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        if (rhocore_f[idx] < 0.0)
            rhocore_f[idx] = 0.0;
        t1 += rhocore_f[idx];
    }


    /* Wait until everyone gets here */
    /*my_barrier(); */

    /*   add the saw-tooth potential induced by applied electric field  */
    init_efield (vnuc_f);


}                               /* end init_nuc */

/******/
