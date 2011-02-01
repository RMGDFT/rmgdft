/************************** SVN Revision Information **************************
 **    $Id: partial_nlop_p.c 779 2007-05-14 19:30:52Z luw $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"

void partial_nlop_p(ION * iptr, REAL * betax, REAL * betay, REAL * betaz,
                    int ip, fftwnd_plan p1, fftwnd_plan p2)
{

    int idx, ix, iy, iz, size;
    REAL ax[3], bx[3], xc, yc, zc;
    REAL r, r3, rsq, rsqd, x, y, z;
    REAL *betax1, *betax2, *betax3, *betay1, *betay2, *betay3, *betaz1, *betaz2, *betaz3;
    REAL hxx, hyy, hzz, cc, t1, t2, invdr;
    fftw_complex *beptr1x, *beptr2x, *beptr3x, *gbptr, *weptr, *gwptr;
    fftw_complex *beptr1y, *beptr2y, *beptr3y, *beptr1z, *beptr2z, *beptr3z;
    SPECIES *sp;

    /* Get species type */
    sp = &ct.sp[iptr->species];
    invdr = ONE / sp->drnlig;

    size = sp->nldim * sp->nldim * sp->nldim;
    my_calloc( weptr, size, fftw_complex );
    my_calloc( gwptr, size, fftw_complex );

    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;
    my_calloc( beptr1x, 3 * size, fftw_complex );
    beptr1y = beptr1x + size;
    beptr1z = beptr1y + size;

    my_calloc( beptr2x, 3 * size, fftw_complex );
    beptr2y = beptr2x + size;
    beptr2z = beptr2y + size;

    my_calloc( beptr3x, 3 * size, fftw_complex );
    beptr3y = beptr3x + size;
    beptr3z = beptr3y + size;
    my_calloc( gbptr, size, fftw_complex );
    if (NULL == gbptr)
        error_handler("can't allocate memerry\n");

    hxx = ct.hxgrid / (REAL) BETA_NX;
    hyy = ct.hygrid / (REAL) BETA_NY;
    hzz = ct.hzgrid / (REAL) BETA_NZ;

    betax1 = betax;
    betax2 = betax1 + ct.max_nlpoints;
    betax3 = betax2 + ct.max_nlpoints;

    betay1 = betay;
    betay2 = betay1 + ct.max_nlpoints;
    betay3 = betay2 + ct.max_nlpoints;

    betaz1 = betaz;
    betaz2 = betaz1 + ct.max_nlpoints;
    betaz3 = betaz2 + ct.max_nlpoints;

    cc = sqrt(3.0 / (4.0 * PI));

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
                to_cartesian(ax, bx);
                t1 = linint(&sp->drbetalig[ip][0], r, invdr);
                t2 = linint(&sp->betalig[ip][0], r, invdr);
                x = bx[0];
                y = bx[1];
                z = bx[2];

                rsq = x * x + y * y + z * z;
                r3 = r * rsq + 1.0e-30;
                rsqd = rsq + 1.0e-20;

                beptr1x[idx].re = cc * ((t2 * (rsq - x * x) / r3) + (x * x * t1 / rsqd));
                beptr1y[idx].re = cc * ((-t2 * x * y / r3) + (x * y * t1 / rsqd));
                beptr1z[idx].re = cc * ((-t2 * x * z / r3) + (x * z * t1 / rsqd));

                beptr2x[idx].re = cc * ((-t2 * x * z / r3) + (x * z * t1 / rsqd));
                beptr2y[idx].re = cc * ((-t2 * y * z / r3) + (y * z * t1 / rsqd));
                beptr2z[idx].re = cc * ((t2 * (rsq - z * z) / r3) + (z * z * t1 / rsqd));

                beptr3x[idx].re = cc * ((-t2 * x * y / r3) + (x * y * t1 / rsqd));
                beptr3y[idx].re = cc * ((t2 * (rsq - y * y) / r3) + (y * y * t1 / rsqd));
                beptr3z[idx].re = cc * ((-t2 * y * z / r3) + (y * z * t1 / rsqd));

                beptr1x[idx].im = 0.0;
                beptr1y[idx].im = 0.0;
                beptr1z[idx].im = 0.0;

                beptr2x[idx].im = 0.0;
                beptr2y[idx].im = 0.0;
                beptr2z[idx].im = 0.0;

                beptr3x[idx].im = 0.0;
                beptr3y[idx].im = 0.0;
                beptr3z[idx].im = 0.0;

                idx++;

                zc += hzz;
            }                   /* end for */

            yc += hyy;
        }                       /* end for */

        xc += hxx;
    }                           /* end for */

    fftwnd_one(p1, beptr1x, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betax1);

    fftwnd_one(p1, beptr2z, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betaz2);

    fftwnd_one(p1, beptr3y, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betay3);

    fftwnd_one(p1, beptr1y, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betay1);
    assign_weight(sp, weptr, betax3);

    fftwnd_one(p1, beptr1z, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betaz1);
    assign_weight(sp, weptr, betax2);

    fftwnd_one(p1, beptr2y, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betay2);
    assign_weight(sp, weptr, betaz3);

    my_free(beptr1x);
    my_free(beptr2x);
    my_free(beptr3x);
    my_free(gbptr);
    my_free(weptr);
    my_free(gwptr);

}
