/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"

void partial_nlop_s(ION * iptr, REAL * betax, REAL * betay, REAL * betaz,
                    int ip, fftwnd_plan p1, fftwnd_plan p2)
{

    int idx, ix, iy, iz, size;
    REAL r, ax[3], bx[3], xc, yc, zc;
    REAL invdr, t1, hxx, hyy, hzz;
    fftw_complex *weptr, *beptrx, *beptry, *beptrz, *gwptr, *gbptr;
    SPECIES *sp;

    /* Get species type */
    sp = &ct.sp[iptr->species];

    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;
    my_calloc( beptrx, 3 * size, fftw_complex );
    beptry = beptrx + size;
    beptrz = beptry + size;
    my_calloc( gbptr, size, fftw_complex );

    size = sp->nldim * sp->nldim * sp->nldim;
    my_calloc( weptr, size, fftw_complex );
    my_calloc( gwptr, size, fftw_complex );
    if (NULL == gwptr)
        error_handler("can't allocate memerry\n");

    hxx = ct.hxgrid / (REAL) BETA_NX;
    hyy = ct.hygrid / (REAL) BETA_NY;
    hzz = ct.hzgrid / (REAL) BETA_NZ;

    invdr = ONE / sp->drnlig;

    idx = 0;
    xc = iptr->nlxcstart;
    for (ix = 0; ix < sp->nlfdim; ix++)
    {

        yc = iptr->nlycstart;
        for (iy = 0; iy < sp->nlfdim; iy++)
        {

            zc = iptr->nlzcstart;
            for (iz = 0; iz < sp->nlfdim; iz++)
            {

                ax[0] = xc - iptr->xtal[0];
                ax[1] = yc - iptr->xtal[1];
                ax[2] = zc - iptr->xtal[2];

                r = metric(ax);
                t1 = linint(&sp->drbetalig[ip][0], r, invdr);
                to_cartesian(ax, bx);
                r += 1.0e-10;

                beptrx[idx].re = sqrt(1.0 / (4.0 * PI)) * t1 * bx[0] / r;
                beptry[idx].re = sqrt(1.0 / (4.0 * PI)) * t1 * bx[1] / r;
                beptrz[idx].re = sqrt(1.0 / (4.0 * PI)) * t1 * bx[2] / r;
                beptrx[idx].im = 0.0;
                beptry[idx].im = 0.0;
                beptrz[idx].im = 0.0;

                idx++;
                zc += hzz;
            }                   /* end for */

            yc += hyy;
        }                       /* end for */

        xc += hxx;
    }                           /* end for */

    fftwnd_one(p1, beptrx, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betax);

    fftwnd_one(p1, beptry, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betay);

    fftwnd_one(p1, beptrz, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betaz);

    my_free(weptr);
    my_free(gwptr);
    my_free(beptrx);
    my_free(gbptr);
}
