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
#include "AtomicInterpolate.h"

void init_nuc (double * vnuc_f, double * rhoc_f, double * rhocore_f)
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

    double r, Zv, rc, rc2, rcnorm, t1;
    double x[3];
    double hxxgrid, hyygrid, hzzgrid;
    double xside, yside, zside;
    SPECIES *sp;
    ION *iptr;
    int npes = get_PE_X() * get_PE_Y() * get_PE_Z();

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
        rhoc_f[idx] = ct.background_charge / FP0_BASIS / get_vel_f() / npes;
        rhocore_f[idx] = 0.0;
        vnuc_f[idx] = 0.0;
    }

    /* Loop over ions determine num of ions whose local pp has overlap with this pe */

    pct.num_loc_ions = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];

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


        bool map_x = false, map_y = false, map_z = false;
        for (ix = xstart; ix < xend; ix++)
        {
            // fold the grid into the unit cell
            // maxium fold 20 times, local potential can extend to 20-1 unit cell
            ixx = (ix + 20 * FNX_GRID) % FNX_GRID;
            if(ixx >= ilow && ixx < ihi)
            {
                map_x = true;
                break;
            }
        }

        if(!map_x) continue;

        for (iy = ystart; iy < yend; iy++)
        {
            // fold the grid into the unit cell
            iyy = (iy + 20 * FNY_GRID) % FNY_GRID;
            if(iyy >= jlow && iyy < jhi)
            {
                map_y = true;
                break;
            }
        }

        if(!map_y) continue;


        for (iz = zstart; iz < zend; iz++)
        {
            // fold the grid into the unit cell
            izz = (iz + 20 * FNZ_GRID) % FNZ_GRID;
            if(izz >= klow && izz < khi)
            {
                map_z = true;
                break;
            }
        }

        if(!map_z) continue;

        pct.loc_ions_list[pct.num_loc_ions] = ion;
        pct.num_loc_ions ++;
    }

    if(pct.localpp != NULL) 
    {
        free(pct.localpp);
        free(pct.localrhoc);
        free(pct.localrhonlcc);
    }

    pct.localpp = (double *) malloc(pct.num_loc_ions * FP0_BASIS * sizeof(double) + 1024);
    pct.localrhoc = (double *) malloc(pct.num_loc_ions * FP0_BASIS * sizeof(double) + 1024);
    pct.localrhonlcc = (double *) malloc(pct.num_loc_ions * FP0_BASIS * sizeof(double) + 1024);
    for(idx = 0; idx < pct.num_loc_ions *FP0_BASIS; idx++)
    {
        pct.localpp[idx] = 0.0;
        pct.localrhoc[idx] = 0.0;
        pct.localrhonlcc[idx] = 0.0;
    }


    int ion1;
    for (ion1 = 0; ion1 < pct.num_loc_ions; ion1++)
    {
        ion = pct.loc_ions_list[ion1];
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


        for (ix = xstart; ix < xend; ix++)
        {
            // fold the grid into the unit cell
            // maxium fold 20 times, local potential can extend to 20-1 unit cell
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

                                t1= AtomicInterpolateInline (&sp->localig[0], r);
                                vnuc_f[idx] += t1;
                                pct.localpp[ion1 * FP0_BASIS + idx] += t1;

                                t1= Zv * exp (-r * r / rc2) * rcnorm;
                                rhoc_f[idx] += t1;
                                pct.localrhoc[ion1 * FP0_BASIS + idx] += t1;
                                

                                if (sp->nlccflag)
                                {
                    
                                    t1 = AtomicInterpolateInline (&sp->rhocorelig[0], r);
                                    rhocore_f[idx] += t1;
                                    pct.localrhonlcc[ion1 * FP0_BASIS + idx] += t1;

                                }

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
    if(ct.runflag == 5 || ct.runflag == 6 || ct.forceflag == TDDFT) return;
    init_efield (vnuc_f);


}                               /* end init_nuc */

/******/
