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
 *   void nlccforce(rmg_double_t *rho, rmg_double_t *vxc)
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



void nlccforce (rmg_double_t * rho, rmg_double_t * vxc)
{

    int ix, iy, iz, ion, idx;
    int *pvec, docount, ishift;
    int ilow, jlow, klow, ihi, jhi, khi, map;
    int *Aix, *Aiy, *Aiz;
    rmg_double_t r, xc, yc, zc, invdr;
    rmg_double_t ax[3], axs[3], bx[3];
    rmg_double_t shift[4];
    rmg_double_t fx, fy, fz;
    SPECIES *sp;
    ION *iptr;
    rmg_double_t deltac;
    rmg_double_t *locsum, *rx, *ry, *rz, *prjptr, *pptr;
    rmg_double_t sumxc2, sumxyc;
#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif


    my_calloc( Aix, FNX_GRID, int );
    my_calloc( Aiy, FNY_GRID, int );
    my_calloc( Aiz, FNZ_GRID, int );

    my_malloc (locsum, 12, rmg_double_t);
    my_malloc (prjptr, 12 * get_FP0_BASIS(), rmg_double_t);
    my_malloc (pvec, get_FP0_BASIS(), int);


    rx = prjptr;
    ry = rx + 4 * get_FP0_BASIS();
    rz = ry + 4 * get_FP0_BASIS();



    deltac = ct.hmaxgrid / 200.0 / (rmg_double_t) FG_NX;
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

            invdr = ONE / sp->drlig;

            /* Determine mapping indices or even if a mapping exists */
            map = get_index (pct.gridpe, iptr, Aix, Aiy, Aiz, &ilow, &ihi, &jlow, &jhi, &klow, &khi,
                             sp->ldim, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(),
                             ct.psi_fnxgrid, ct.psi_fnygrid, ct.psi_fnzgrid,
                             &iptr->lxcstart, &iptr->lycstart, &iptr->lzcstart);


            /* If there is any overlap then we have to generate the mapping */
            /*my_barrier(); */
            docount = 0;
            if (map)
            {

                zc = iptr->lzcstart;
                for (iz = 0; iz < sp->ldim; iz++)
                {

                    yc = iptr->lycstart;
                    for (iy = 0; iy < sp->ldim; iy++)
                    {

                        xc = iptr->lxcstart;
                        for (ix = 0; ix < sp->ldim; ix++)
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


                                for (ishift = 0; ishift < 4; ishift++)
                                {

                                    axs[0] = bx[0] - shift[ishift];
                                    axs[1] = bx[1];
                                    axs[2] = bx[2];
                                    r = sqrt (axs[0] * axs[0] + axs[1] * axs[1] + axs[2] * axs[2]);

                                    rx[docount + ishift * get_FP0_BASIS()] =
                                        linint (&sp->rhocorelig[0], r, invdr);



                                }       /* end for */


                                for (ishift = 0; ishift < 4; ishift++)
                                {

                                    axs[0] = bx[0];
                                    axs[1] = bx[1] - shift[ishift];
                                    axs[2] = bx[2];
                                    r = sqrt (axs[0] * axs[0] + axs[1] * axs[1] + axs[2] * axs[2]);

                                    ry[docount + ishift * get_FP0_BASIS()] =
                                        linint (&sp->rhocorelig[0], r, invdr);


                                }       /* end for */

                                for (ishift = 0; ishift < 4; ishift++)
                                {

                                    axs[0] = bx[0];
                                    axs[1] = bx[1];
                                    axs[2] = bx[2] - shift[ishift];
                                    r = sqrt (axs[0] * axs[0] + axs[1] * axs[1] + axs[2] * axs[2]);

                                    rz[docount + ishift * get_FP0_BASIS()] =
                                        linint (&sp->rhocorelig[0], r, invdr);



                                }       /* end for */



                                docount++;


                            }   /* end if */

                            xc += ct.hxxgrid;

                        }       /* end for */

                        yc += ct.hyygrid;

                    }           /* end for */

                    zc += ct.hzzgrid;

                }               /* end for */



            }                   /* end if */

            /*my_barrier(); */


            pptr = prjptr;
            for (ishift = 0; ishift < 12; ishift++)
            {

                locsum[ishift] = ZERO;

                for (idx = 0; idx < docount; idx++)
                {

                    locsum[ishift] += pptr[idx] * vxc[pvec[idx]];

                }               /* end for */

                locsum[ishift] = locsum[ishift] * ct.vel_f;
                pptr += get_FP0_BASIS();

            }                   /* end for */


            idx = 12;
            global_sums (locsum, &idx, pct.img_comm);
	    if (ct.spin_flag)
	    	for (ishift = 0; ishift < 12; ishift++)
			locsum[ishift] *= 0.5;         
	        /* factor 0.5 is because when calculating exchange correlation
		half of nonlinear core corection charge is added to spin up and down density */


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
        }                       /* end if */


    }                           /* end for */



    my_free (Aiz);
    my_free (Aiy);
    my_free (Aix);

    my_free (pvec);
    my_free (prjptr);
    my_free (locsum);

#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (NLCCFORCE_TIME, (time2 - time1));
#endif

}                               /* end nlccforce */

/******/
