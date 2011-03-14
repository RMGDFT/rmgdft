/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/lforce.c *****
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"



void lforce(REAL * rho, REAL * vh)
{

    int ix, iy, iz;
    int ion, idx;
    int *pvec, icount;
    int ilow, jlow, klow, ihi, jhi, khi, map;
    int *Aix, *Aiy, *Aiz;
    REAL r, xc, yc, zc, Zv, rc, rcnorm, t1;
    REAL x, y, z, rc2, invdr, norm1;
    REAL fx, fy, fz;
    SPECIES *sp;
    ION *iptr;
    REAL *rx, *ry, *rz;
    REAL *urx, *ury, *urz;
    int L0_LDIM, ixstart, iystart, izstart, ii, jj, kk;
    REAL hxgrid, hygrid, hzgrid;

    REAL time1, time2;
    time1 = my_crtc();

    /* Grab some memory for temporary storage */
    /*shuchun wang */
    my_calloc( pvec, FP0_BASIS, int );
    my_calloc( Aix, FNX_GRID, int );
    my_calloc( Aiy, FNY_GRID, int );
    my_calloc( Aiz, FNZ_GRID, int );


    my_malloc_init( rx, 6 * FP0_BASIS, REAL );
    ry = rx + FP0_BASIS;
    rz = ry + FP0_BASIS;
    urx = rz + FP0_BASIS;
    ury = urx + FP0_BASIS;
    urz = ury + FP0_BASIS;


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
        rcnorm = rc * rc * rc * pow(PI, 1.5);
        rcnorm = 1.0 / rcnorm;
        norm1 = -Zv * rcnorm / rc2;

        L0_LDIM = sp->ldim;

        hxgrid = ct.hxxgrid * ct.xside;
        hygrid = ct.hyygrid * ct.yside;
        hzgrid = ct.hzzgrid * ct.zside;

        /* Generate range of indices over which the short-range difference */
        /* potential will be mapped onto the global grid.                  */
        get_start(L0_LDIM, ct.ions[ion].crds[0], ct.xcstart, hxgrid, &ixstart, &iptr->lxcstart);
        get_start(L0_LDIM, ct.ions[ion].crds[1], ct.ycstart, hygrid, &iystart, &iptr->lycstart);
        get_start(L0_LDIM, ct.ions[ion].crds[2], ct.zcstart, hzgrid, &izstart, &iptr->lzcstart);

        /* Next we have to map these onto the global grid */

        /* Generate indices */
        get_Ai(L0_LDIM, Aix, FNX_GRID, ixstart);
        get_Ai(L0_LDIM, Aiy, FNY_GRID, iystart);
        get_Ai(L0_LDIM, Aiz, FNZ_GRID, izstart);


        /* Now we need to determine if any of this ions local short */
        /* ranged difference potential maps onto this particular    */
        /* processors space.                                        */
        pe2xyz(pct.gridpe, &ii, &jj, &kk);
        ilow = ii * FPX0_GRID;
        jlow = jj * FPY0_GRID;
        klow = kk * FPZ0_GRID;
        ihi = ilow + FPX0_GRID - 1;
        jhi = jlow + FPY0_GRID - 1;
        khi = klow + FPZ0_GRID - 1;

        ii = jj = kk = FALSE;
        for (idx = 0; idx < L0_LDIM; idx++)
        {

            if ((Aix[idx] >= ilow) && (Aix[idx] <= ihi))
                ii = TRUE;
            if ((Aiy[idx] >= jlow) && (Aiy[idx] <= jhi))
                jj = TRUE;
            if ((Aiz[idx] >= klow) && (Aiz[idx] <= khi))
                kk = TRUE;

        }                       /* end for */

        map = ii & jj;
        map = map & kk;

        icount = 0;
        /* If there is any overlap then we have to generate the mapping */
        if (map)
        {

            invdr = 1.0 / sp->drlig;

            xc = iptr->lxcstart;
            for (ix = 0; ix < L0_LDIM; ix++)
            {

                yc = iptr->lycstart;
                for (iy = 0; iy < L0_LDIM; iy++)
                {

                    zc = iptr->lzcstart;
                    for (iz = 0; iz < L0_LDIM; iz++)
                    {


                        if (((Aix[ix] >= ilow) && (Aix[ix] <= ihi)) &&
                            ((Aiy[iy] >= jlow) && (Aiy[iy] <= jhi)) &&
                            ((Aiz[iz] >= klow) && (Aiz[iz] <= khi)))
                        {

                            pvec[icount] =
                                FPY0_GRID * FPZ0_GRID * (Aix[ix] % FPX0_GRID) +
                                FPZ0_GRID * (Aiy[iy] % FPY0_GRID) + (Aiz[iz] % FPZ0_GRID);

                            x = xc - iptr->crds[0];
                            y = yc - iptr->crds[1];
                            z = zc - iptr->crds[2];

                            r = sqrt(x * x + y * y + z * z);
                            t1 = 2.0 * norm1 * exp(-r * r / rc2);
                            urx[icount] = x * t1;
                            ury[icount] = y * t1;
                            urz[icount] = z * t1;

                            t1 = linint(sp->drlocalig, r, invdr);
                            r += 1.0e-10;

                            rx[icount] = -t1 * x / r;
                            ry[icount] = -t1 * y / r;
                            rz[icount] = -t1 * z / r;

                            icount++;

                        }

                        zc += hzgrid;

                    }           /* end for iz */

                    yc += hygrid;

                }               /* end for iy */

                xc += hxgrid;

            }                   /* end for ix */


        }                       /*end if(map) */

        my_barrier();

        fx = fy = fz = 0.0;
        for (idx = 0; idx < icount; idx++)
        {

            fx += rx[idx] * rho[pvec[idx]];
            fy += ry[idx] * rho[pvec[idx]];
            fz += rz[idx] * rho[pvec[idx]];
            fx += urx[idx] * vh[pvec[idx]];
            fy += ury[idx] * vh[pvec[idx]];
            fz += urz[idx] * vh[pvec[idx]];

        }

        iptr->force[ct.fpt[0]][0] -= ct.vel_f * real_sum_all(fx);
        iptr->force[ct.fpt[0]][1] -= ct.vel_f * real_sum_all(fy);
        iptr->force[ct.fpt[0]][2] -= ct.vel_f * real_sum_all(fz);

    }                           /* end for ion */

    /* Release our memory */
    my_free(rx);

    my_free(pvec);
    my_free(Aix);
    my_free(Aiy);
    my_free(Aiz);

    time2 = my_crtc();
    rmg_timings(LFORCE_TIME, time2 - time1, 0);

}

/******/
