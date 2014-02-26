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



void lforce (rmg_double_t * rho, rmg_double_t * vh)
{

    int ix, iy, iz;
    int ion, idx;
    int *pvec, docount, size;
    int ilow, jlow, klow, ihi, jhi, khi, map;
    int icut, icenter;
    int *Aix, *Aiy, *Aiz;
    rmg_double_t r, xc, yc, zc, Zv, rc, rcnorm, t1;
    rmg_double_t ax[3], bx[3], x_hat, y_hat, z_hat, rc2, invdr, norm1;
    rmg_double_t hxxgrid, hyygrid, hzzgrid;
    rmg_double_t fx, fy, fz;
    SPECIES *sp;
    ION *iptr;
    rmg_double_t *rx, *ry, *rz;
    rmg_double_t *urx, *ury, *urz;
#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif

    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();

    my_calloc( Aix, get_FNX_GRID(), int );
    my_calloc( Aiy, get_FNY_GRID(), int );
    my_calloc( Aiz, get_FNZ_GRID(), int );

    size = get_FP0_BASIS();
    my_malloc (rx, 6 * size, rmg_double_t);
    ry = rx + size;
    rz = ry + size;
    urx = rz + size;
    ury = urx + size;
    urz = ury + size;
    my_malloc (pvec, size, int);


    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Generate ion pointer */
        iptr = &ct.ions[ion];


        /* Get species type */
        sp = &ct.sp[iptr->species];

        icenter = sp->ldim / 2;
        icut = (icenter + 1) * (icenter + 1);


        Zv = sp->zvalence;
        rc = sp->rc;
        rc2 = rc * rc;
        rcnorm = rc * rc * rc * pow (PI, 1.5);
        rcnorm = ONE / rcnorm;
        rc2 = sp->rc * sp->rc;
        norm1 = -Zv * rcnorm / rc2;
        invdr = ONE / sp->drlig;


        /* Determine mapping indices or even if a mapping exists */
        map = get_index (pct.gridpe, iptr, Aix, Aiy, Aiz, &ilow, &ihi, &jlow, &jhi, &klow, &khi,
                         sp->ldim, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(),
                         get_FNX_GRID(), get_FNY_GRID(), get_FNZ_GRID(),
                         &iptr->lxcstart, &iptr->lycstart, &iptr->lzcstart);



        /* If there is any overlap then we have to generate the mapping */
        /*my_barrier(); */
        docount = 0;
        if (map)
        {

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


                            pvec[docount] =
                                get_FPY0_GRID() * get_FPZ0_GRID() * ((Aix[ix]-get_FPX_OFFSET()) % get_FPX0_GRID()) +
                                get_FPZ0_GRID() * ((Aiy[iy]-get_FPY_OFFSET()) % get_FPY0_GRID()) +
                                ((Aiz[iz]-get_FPZ_OFFSET()) % get_FPZ0_GRID());


                            ax[0] = xc - iptr->xtal[0];
                            ax[1] = yc - iptr->xtal[1];
                            ax[2] = zc - iptr->xtal[2];
                            to_cartesian (ax, bx);
                            r = metric (ax);

                            t1 = 2.0 * norm1 * exp (-r * r / rc2);
                            urx[docount] = bx[0] * t1;
                            ury[docount] = bx[1] * t1;
                            urz[docount] = bx[2] * t1;

                            t1 = linint (&sp->drlocalig[0], r, invdr);
                            r += 1.0e-10;
                            x_hat = bx[0] / r;
                            y_hat = bx[1] / r;
                            z_hat = bx[2] / r;


                            rx[docount] = -t1 * x_hat;
                            ry[docount] = -t1 * y_hat;
                            rz[docount] = -t1 * z_hat;

                            docount++;


                        }       /* end if */

                        zc += hzzgrid;

                    }           /* end for */

                    yc += hyygrid;

                }               /* end for */

                xc += hxxgrid;

            }                   /* end for */

        }                       /* end if */

        /*my_barrier(); */


        fx = fy = fz = 0.0;
        for (idx = 0; idx < docount; idx++)
        {

            fx += rx[idx] * rho[pvec[idx]];
            fy += ry[idx] * rho[pvec[idx]];
            fz += rz[idx] * rho[pvec[idx]];
            fx += urx[idx] * vh[pvec[idx]];
            fy += ury[idx] * vh[pvec[idx]];
            fz += urz[idx] * vh[pvec[idx]];

        }                       /* end for */

        iptr->force[ct.fpt[0]][0] -= get_vel_f() * real_sum_all (fx, pct.grid_comm);
        iptr->force[ct.fpt[0]][1] -= get_vel_f() * real_sum_all (fy, pct.grid_comm);
        iptr->force[ct.fpt[0]][2] -= get_vel_f() * real_sum_all (fz, pct.grid_comm);

    }                           /* end for */


    my_free(Aiz);
    my_free(Aiy);
    my_free(Aix);
    my_free (pvec);
    my_free (rx);


#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (LFORCE_TIME, (time2 - time1));
#endif

}                               /* end lforce */

/******/
