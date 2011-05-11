/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/nlccforce.c *****/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"



void nlccforce(REAL * rho, REAL * vxc)
{

    int ix, iy, iz, ion, idx;
    int *pvec, icount, ishift;
    int ilow, jlow, klow, ihi, jhi, khi, map;
    int *Aix, *Aiy, *Aiz;
    REAL r, xc, yc, zc, invdr;
    REAL axs[3], ax[3];
    REAL shift[4];
    REAL fx, fy, fz;
    SPECIES *sp;
    ION *iptr;
    REAL deltac;
    REAL *locsum, *rx, *ry, *rz, *prjptr, *pptr;
    REAL sumxc2, sumxyc;
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

    my_malloc_init( locsum, 12, REAL );
    my_malloc_init( prjptr, 12 * FP0_BASIS, REAL );
    rx = prjptr;
    ry = rx + 4 * FP0_BASIS;
    rz = ry + 4 * FP0_BASIS;



    deltac = ct.hmaxgrid / 200.0 / (REAL) RHO_NX;
    shift[0] = -2.0 * deltac;
    shift[1] = 2.0 * deltac;
    shift[2] = -deltac;
    shift[3] = deltac;

    sumxc2 = 2.0 * (shift[0] * shift[0] + shift[2] * shift[2]);

    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];

        if (sp->nlccflag)
        {

            invdr = 1.0 / sp->drlig;

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

            }                   /* end for */

            map = ii & jj;
            map = map & kk;

            icount = 0;
            /* If there is any overlap then we have to generate the mapping */
            if (map)
            {

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
                                    FPY0_GRID * FPZ0_GRID * (Aix[ix] % FPX0_GRID) +
                                    FPZ0_GRID * (Aiy[iy] % FPY0_GRID) + (Aiz[iz] % FPZ0_GRID);

                                ax[0] = xc - iptr->crds[0];
                                ax[1] = yc - iptr->crds[1];
                                ax[2] = zc - iptr->crds[2];


                                for (ishift = 0; ishift < 4; ishift++)
                                {

                                    axs[0] = ax[0] - shift[ishift];
                                    axs[1] = ax[1];
                                    axs[2] = ax[2];
                                    r = sqrt(axs[0] * axs[0] + axs[1] * axs[1] + axs[2] * axs[2]);

                                    rx[icount + ishift * FP0_BASIS] =
                                        linint(&sp->rhocorelig[0], r, invdr);



                                }       /* end for */


                                for (ishift = 0; ishift < 4; ishift++)
                                {

                                    axs[0] = ax[0];
                                    axs[1] = ax[1] - shift[ishift];
                                    axs[2] = ax[2];
                                    r = sqrt(axs[0] * axs[0] + axs[1] * axs[1] + axs[2] * axs[2]);

                                    ry[icount + ishift * FP0_BASIS] =
                                        linint(&sp->rhocorelig[0], r, invdr);


                                }       /* end for */

                                for (ishift = 0; ishift < 4; ishift++)
                                {

                                    axs[0] = ax[0];
                                    axs[1] = ax[1];
                                    axs[2] = ax[2] - shift[ishift];
                                    r = sqrt(axs[0] * axs[0] + axs[1] * axs[1] + axs[2] * axs[2]);

                                    rz[icount + ishift * FP0_BASIS] =
                                        linint(&sp->rhocorelig[0], r, invdr);



                                }



                                icount++;


                            }

                            xc += hxgrid;

                        }       /* end for ix */

                        yc += hygrid;

                    }           /* end for iy */

                    zc += hzgrid;

                }               /* end for iz */



            }                   /* end if map */

            my_barrier();


            pptr = prjptr;
            for (ishift = 0; ishift < 12; ishift++)
            {

                locsum[ishift] = ZERO;

                for (idx = 0; idx < icount; idx++)
                {

                    locsum[ishift] += pptr[idx] * vxc[pvec[idx]];

                }               /* end for */

                locsum[ishift] = locsum[ishift] * ct.vel_f;
                pptr += FP0_BASIS;

            }                   /* end for */


            idx = 12;
            global_sums(locsum, &idx);


            sumxyc = ZERO;
            sumxyc += locsum[0] * shift[0];
            sumxyc += locsum[1] * shift[1];
            sumxyc += locsum[2] * shift[2];
            sumxyc += locsum[3] * shift[3];
            fx = -sumxyc / sumxc2;
            iptr->force[ct.fpt[0]][0] += fx;

            sumxyc = ZERO;
            sumxyc += locsum[4] * shift[0];
            sumxyc += locsum[5] * shift[1];
            sumxyc += locsum[6] * shift[2];
            sumxyc += locsum[7] * shift[3];
            fy = -sumxyc / sumxc2;
            iptr->force[ct.fpt[0]][1] += fy;

            sumxyc = ZERO;
            sumxyc += locsum[8] * shift[0];
            sumxyc += locsum[9] * shift[1];
            sumxyc += locsum[10] * shift[2];
            sumxyc += locsum[11] * shift[3];
            fz = -sumxyc / sumxc2;
            iptr->force[ct.fpt[0]][2] += fz;

        }                       /* end if sp->nlccflag */


    }                           /* end for ion */



    my_free(prjptr);
    my_free(locsum);

    my_free(pvec);
    my_free(Aix);
    my_free(Aiy);
    my_free(Aiz);

    time2 = my_crtc();
    rmg_timings(NLCCFORCE_TIME, time2 - time1);

}

/******/
