/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*


                            init_nuc.c


    Set up the nuclear local potential on the finest grid and the gaussian
    compensating charges.



*/



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"

#define DELTA		1e-10
#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.12837916709551257390  /* 2./sqrt(pi) */
#endif
double erf(double);


static void init_vcomp(double *);
static double get_charge(double *);




void init_nuc(double *vnuc, double *rhoc, double *rhocore)
{

    int ix, iy, iz;
    int ion, idx;
    int ii, jj, kk;
    int ilow, jlow, klow, ihi, jhi, khi, map;
    int ixstart, iystart, izstart;
    int *Aix, *Aiy, *Aiz;
    int icount;
    int *pvec;
    int number = FP0_BASIS;
    int storage_spacing = 1;
    double r, r2, xc, yc, zc, Zv, rc, rcnorm, t1;
    double x, y, z, rc2, invdr, charge_bg, x_locate;
    SPECIES *sp;
    ION *iptr;

    int L0_LDIM;
    REAL hxgrid, hygrid, hzgrid;

    if (pct.gridpe == 0)
    {

        printf(" Begin init_nuc ...\n");

    }                           /* end if */

    fflush(NULL);

    my_barrier();


    /* Grab some memory for temporary storage */
    /*shuchun wang */
    my_calloc( pvec, FP0_BASIS, int );
    my_calloc( Aix, FNX_GRID, int );
    my_calloc( Aiy, FNY_GRID, int );
    my_calloc( Aiz, FNZ_GRID, int );

    /* Zero out the nuclear potential array */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        vnuc[idx] = 0.;

    }                           /* end for */

        pe2xyz(pct.gridpe, &ii, &jj, &kk);
        ilow = pct.FPX_OFFSET;
        jlow = pct.FPY_OFFSET;
        klow = pct.FPZ_OFFSET;
        ihi = ilow + FPX0_GRID - 1;
        jhi = jlow + FPY0_GRID - 1;
        khi = klow + FPZ0_GRID - 1;


    /* Initialize the compensating charge array and the core charge
       array */
    t1 = ct.background_charge / (double) FP0_BASIS / ct.vel_f / NPES;
    
    printf("\n bg_begin = %f and bg_end = %f  BT = %f  t1=%f \n", ct.bg_begin, ct.bg_end, ct.BT, t1);   


    if (ct.background_charge == 0)
    {
	    for (idx = 0; idx < FP0_BASIS; idx++)
	    {

		    rhoc[idx] = t1;

	    } 
    }
    else
    {
	    for (ix = 0; ix < FPX0_GRID; ix++)
	    {

		    x_locate = (ix + ilow)/2.0; // everything is based on coarse grid now!
		    charge_bg = t1 / (1.0 + exp((ct.bg_begin - x_locate)/ct.BT)  + exp((x_locate - ct.bg_end)/ct.BT) ); 


		    for (iy = 0; iy < FPY0_GRID; iy++)
		    {

			    for (iz = 0; iz < FPZ0_GRID; iz++)
			    {

				    rhoc[ix * FPY0_GRID * FPZ0_GRID + iy * FPZ0_GRID + iz] = charge_bg;


			    }                   /* end for */

		    }                       /* end for */


		    if (jj == 0 && kk == 0)
			    printf("x_locate = %5f charge_bg  = %15.12f \n ", x_locate, charge_bg );

	    }                           /* end for */

	    /* Check background charges */
	    ct.crho = get_charge(rhoc);



	    if (pct.gridpe == 0)
	    {
		    printf("\n background charges: %e \n", ct.crho);
		    printf(" Rescaling background charges\n");
	    }
	    t1 = ct.background_charge / ct.crho;
	    sscal(&number, &t1, rhoc, &storage_spacing);

	    ct.crho = get_charge(rhoc);

	    if (pct.gridpe == 0)
		    printf(" Rescaled compensating charges: %e \n", ct.crho);
    }




    for (idx = 0; idx < FP0_BASIS; idx++)
    {

	    rhocore[idx] = 0.;

    }                           /* end for */



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
        rcnorm = 1. / rcnorm;

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


        ii = jj = kk = FALSE;
        /* Now we need to determine if any of this ions local short */
        /* ranged difference potential maps onto this particular    */
        /* processors space.                                        */
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

        /* If there is any overlap then we have to generate the mapping */
        if (map)
        {


            invdr = 1. / sp->drlig;
            icount = 0;
            zc = iptr->lzcstart;

            for (iz = 0; iz < L0_LDIM; iz++)
            {

                yc = iptr->lycstart;

                for (iy = 0; iy < L0_LDIM; iy++)
                {

                    xc = iptr->lxcstart;

                    for (ix = 0; ix < L0_LDIM; ix++)
                    {

                        if (((Aix[ix] >= ilow) && (Aix[ix] <= ihi)) &&
                            ((Aiy[iy] >= jlow) && (Aiy[iy] <= jhi)) &&
                            ((Aiz[iz] >= klow) && (Aiz[iz] <= khi)))
                        {

                            pvec[icount] =
                                FPY0_GRID * FPZ0_GRID * (Aix[ix] - pct.FPX_OFFSET) +
                                FPZ0_GRID * (Aiy[iy] - pct.FPY_OFFSET) + (Aiz[iz] - pct.FPZ_OFFSET);

                            x = xc - iptr->crds[0];
                            y = yc - iptr->crds[1];
                            z = zc - iptr->crds[2];

                            r2 = x * x + y * y + z * z;
                            r = sqrt(r2);


                            rhoc[pvec[icount]] += Zv * exp(-r2 / rc2) * rcnorm;

                            vnuc[pvec[icount]] += linint(sp->localig, r, invdr);

                            if (sp->nlccflag)
                                rhocore[pvec[icount]] += linint(sp->rhocorelig, r, invdr);



                            icount++;

                        }       /* end if */

                        xc += hxgrid;

                    }           /* end for */

                    yc += hygrid;

                }               /* end for */

                zc += hzgrid;

            }                   /* end for */


        }                       /* end if */

    }                           /* end for */



    /* Check compensating charges */
    ct.crho = get_charge(rhoc);



    /* Rescale compensating charges */
    if (pct.gridpe == 0)
    {
        printf("\n compensating charges(background_charge + ionic charge): %e \n", ct.crho);
        printf(" Rescaling compensating charges\n");
    }
    t1 = ct.nel / ct.crho;
    if (pct.gridpe == 0)
        printf(" Rescaled compensating charges: ratio -1 = %e \n", t1-1.0);
    sscal(&number, &t1, rhoc, &storage_spacing);

    /* Check new compensating charges */
    ct.crho = get_charge(rhoc);

    if (pct.gridpe == 0)
        printf(" Rescaled compensating charges: %e \n", ct.crho);


    t1 = get_charge(rhocore);

    if (pct.gridpe == 0)
        printf(" core charges: %e \n", t1);


    /* Release our memory */
    my_free(pvec);
    my_free(Aix);
    my_free(Aiy);
    my_free(Aiz);

    if (pct.gridpe == 0)
    {

        printf(" init_nuc done\n");

    }                           /* end if */

    /* Wait until everyone gets here */
    fflush(NULL);
    my_barrier();

}                               /* end init_nuc */


/*
    Initialization of the compensating potential
*/
static void init_vcomp(double *vc)
{
    int ix, iy, iz, ii, jj, kk;
    int ilow, jlow, klow;
    int idx, ion, istart, jstart;
    double r, Zv, rc, point[3];
    SPECIES *sp;
    ION *iptr;


    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        vc[idx] = 0.;

    }                           /* end for */


    ilow = pct.FPX_OFFSET;
    jlow = pct.FPY_OFFSET;
    klow = pct.FPZ_OFFSET;


    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];

        Zv = sp->zvalence;
        rc = sp->rc;

        for (ix = 0; ix < FPX0_GRID; ix++)
        {

            point[0] = (ix + ilow) * ct.hxxgrid;
            istart = FPZ0_GRID * FPY0_GRID * ix;

            for (iy = 0; iy < FPY0_GRID; iy++)
            {

                point[1] = (iy + jlow) * ct.hyygrid;
                jstart = istart + FPZ0_GRID * iy;

                for (iz = 0; iz < FPZ0_GRID; iz++)
                {

                    point[2] = (iz + klow) * ct.hzzgrid;

                    r = minimage1(point, iptr->crds);


                    if (r <= DELTA)
                    {
                        vc[jstart + iz] += Zv * M_2_SQRTPI / rc;
                    }
                    else
                    {
                        vc[jstart + iz] += Zv * erf(r / rc) / r;
                    }

                }               /* end for */

            }                   /* end for */

        }                       /* end for */
    }
}



static double get_charge(double *rho1)
{
    int idx;
    double charge;

    charge = 0.;
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        charge += rho1[idx];

    }                           /* end for */
    charge = real_sum_all(charge, pct.grid_comm);
    charge *= ct.vel_f;

    return charge;
}
