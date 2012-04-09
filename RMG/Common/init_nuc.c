/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"


void init_nuc (REAL * vnuc_f, REAL * rhoc_f, REAL * rhocore_f)
{

    int ix, iy, iz;
    int ion, idx;
    int ilow, jlow, klow, ihi, jhi, khi, map;
    int Aix[FNX_GRID], Aiy[FNY_GRID], Aiz[FNZ_GRID];
    int icount;
    int *pvec;
    REAL r, xc, yc, zc, Zv, rc, rc2, rcnorm, t1;
    REAL x[3], invdr;
    SPECIES *sp;
    ION *iptr;

    /* Grab some memory for temporary storage */
    my_malloc (pvec, FP0_BASIS, int);

    /* Initialize the compensating charge array and the core charge array */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        rhoc_f[idx] = ct.background_charge / FP0_BASIS / ct.vel_f / NPES;
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

        /* Determine mapping indices or even if a mapping exists */
        map = get_index (pct.gridpe, iptr, Aix, Aiy, Aiz, &ilow, &ihi, &jlow, &jhi, &klow, &khi,
                         sp->ldim, FPX0_GRID, FPY0_GRID, FPZ0_GRID,
                         ct.psi_fnxgrid, ct.psi_fnygrid, ct.psi_fnzgrid,
                         &iptr->lxcstart, &iptr->lycstart, &iptr->lzcstart);


	icount = 0;
        /* If there is any overlap then we have to generate the mapping */
        if (map)
        {
            invdr = 1.0 / sp->drlig;
            icount = 0;

            xc = iptr->lxcstart;
            for (ix = 0; ix < sp->ldim; ix++)
            {
                yc = iptr->lycstart;
                for (iy = 0; iy < sp->ldim; iy++)
                {
                    zc = iptr->lzcstart;
                    for (iz = 0; iz < sp->ldim; iz++)
                    {
                        if (((Aix[ix] >= ilow) && (Aix[ix] <= ihi)) &&
                            ((Aiy[iy] >= jlow) && (Aiy[iy] <= jhi)) &&
                            ((Aiz[iz] >= klow) && (Aiz[iz] <= khi)))
                        {
                            pvec[icount] =
                                FPY0_GRID * FPZ0_GRID * (Aix[ix] % FPX0_GRID) +
                                FPZ0_GRID * (Aiy[iy] % FPY0_GRID) + (Aiz[iz] % FPZ0_GRID);

                            x[0] = xc - iptr->xtal[0];
                            x[1] = yc - iptr->xtal[1];
                            x[2] = zc - iptr->xtal[2];
                            r = metric (x);

                            vnuc_f[pvec[icount]] += linint (&sp->localig[0], r, invdr);
                            rhoc_f[pvec[icount]] += Zv * exp (-r * r / rc2) * rcnorm;

                            if (sp->nlccflag)
                                rhocore_f[pvec[icount]] += linint (&sp->rhocorelig[0], r, invdr);


                            icount++;
                        }

                        zc += ct.hzzgrid;

                    }           /* end for */

                    yc += ct.hyygrid;

                }               /* end for */

                xc += ct.hxxgrid;

            }                   /* end for */

        }                       /* end if */

        pct.lptrlen[ion] = icount;

    }                           /* end for */

    /* Check compensating charges */
    ct.crho = 0.0;
    for (idx = 0; idx < FP0_BASIS; idx++)
        ct.crho += rhoc_f[idx];


    ct.crho = ct.crho * ct.vel_f;
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


    /*t1 = t1 * ct.vel_f;
    t1 = real_sum_all (t1, pct.grid_comm);*/


    /* Release our memory */
    my_free (pvec);


    /* Wait until everyone gets here */
    /*my_barrier(); */

    /*   add the saw-tooth potential induced by applied electric field  */
    init_efield (vnuc_f);


}                               /* end init_nuc */

/******/
