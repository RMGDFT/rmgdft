/************************** SVN Revision Information **************************
 **    $Id: partial_nlop_d.c 779 2007-05-14 19:30:52Z luw $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"

void partial_nlop_d(ION * iptr, REAL * betax, REAL * betay, REAL * betaz,
                    int ip, fftwnd_plan p1, fftwnd_plan p2)
{

    int idx, ix, iy, iz, size;
    REAL ax[3], bx[3], xc, yc, zc, t1, t2, t3;
    REAL cc, hxx, hyy, hzz, invdr, x, y, z;
    REAL dy1_dx, dy1_dy, dy1_dz;
    REAL dy2_dx, dy2_dy, dy2_dz;
    REAL dy3_dx, dy3_dy, dy3_dz;
    REAL dy4_dx, dy4_dy, dy4_dz;
    REAL dy5_dx, dy5_dy, dy5_dz;
    REAL dt2_dx, dt2_dy, dt2_dz;
    REAL y1, y2, y3, y4, y5;
    REAL xsq, ysq, zsq, r, rsq, rsqd, r3, r4;
    fftw_complex *beptr1x, *beptr2x, *beptr3x, *beptr4x, *beptr5x, *weptr, *gbptr, *gwptr;
    fftw_complex *beptr1y, *beptr2y, *beptr3y, *beptr4y, *beptr5y;
    fftw_complex *beptr1z, *beptr2z, *beptr3z, *beptr4z, *beptr5z;
    REAL *betax1, *betax2, *betax3, *betax4, *betax5;
    REAL *betay1, *betay2, *betay3, *betay4, *betay5;
    REAL *betaz1, *betaz2, *betaz3, *betaz4, *betaz5;
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

    my_calloc( beptr4x, 3 * size, fftw_complex );
    beptr4y = beptr4x + size;
    beptr4z = beptr4y + size;

    my_calloc( beptr5x, 3 * size, fftw_complex );
    beptr5y = beptr5x + size;
    beptr5z = beptr5y + size;

    my_calloc( gbptr, size, fftw_complex );
    if (NULL == gbptr)
        error_handler("can't allocate memerry\n");


    hxx = ct.hxgrid / (REAL) BETA_NX;
    hyy = ct.hygrid / (REAL) BETA_NY;
    hzz = ct.hzgrid / (REAL) BETA_NZ;

    betax1 = betax;
    betax2 = betax1 + ct.max_nlpoints;
    betax3 = betax2 + ct.max_nlpoints;
    betax4 = betax3 + ct.max_nlpoints;
    betax5 = betax4 + ct.max_nlpoints;

    betay1 = betay;
    betay2 = betay1 + ct.max_nlpoints;
    betay3 = betay2 + ct.max_nlpoints;
    betay4 = betay3 + ct.max_nlpoints;
    betay5 = betay4 + ct.max_nlpoints;

    betaz1 = betaz;
    betaz2 = betaz1 + ct.max_nlpoints;
    betaz3 = betaz2 + ct.max_nlpoints;
    betaz4 = betaz3 + ct.max_nlpoints;
    betaz5 = betaz4 + ct.max_nlpoints;

    cc = sqrt(5.0 / (4.0 * PI));
    t3 = sqrt(3.0);

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
                x = bx[0];
                y = bx[1];
                z = bx[2];

                rsq = x * x + y * y + z * z;
                r3 = r * rsq + 1.0e-30;
                rsqd = rsq + 1.0e-20;
                r4 = rsq * rsq + 1.e-40;

                t1 = t3 * linint(&sp->drbetalig[ip][0], r, invdr);
                t2 = t3 * linint(&sp->betalig[ip][0], r, invdr);
                r += 1.0e-10;


                xsq = x * x;
                ysq = y * y;
                zsq = z * z;

                y1 = x * y / rsqd;
                y2 = z * x / rsqd;
                y3 = (t3 * z * z - rsq / t3) / (2.0 * rsqd);
                y4 = y * z / rsqd;
                y5 = (x * x - y * y) / (2.0 * rsqd);

                dy1_dx = (y - 2.0 * xsq * y / rsqd) / rsqd;
                dy1_dy = (x - 2.0 * x * ysq / rsqd) / rsqd;
                dy1_dz = -2.0 * x * y * z / r4;

                dy2_dx = (z - 2.0 * z * xsq / rsqd) / rsqd;
                dy2_dy = -2.0 * x * y * z / r4;
                dy2_dz = (x - 2.0 * zsq * x / rsqd) / rsqd;

                dy3_dx = -t3 * x * z * z / (rsqd * rsqd);
                dy3_dy = -t3 * y * z * z / (rsqd * rsqd);
                dy3_dz = t3 * (z - z * zsq / rsqd) / rsqd;

                dy4_dx = -2.0 * x * y * z / r4;
                dy4_dy = (z - 2.0 * z * ysq / rsqd) / rsqd;
                dy4_dz = (y - 2.0 * zsq * y / rsqd) / rsqd;

                dy5_dx = x * (1 - (x * x - y * y) / rsqd) / rsqd;
                dy5_dy = y * (-1 - (x * x - y * y) / rsqd) / rsqd;
                dy5_dz = -z * (x * x - y * y) / (rsqd * rsqd);

                dt2_dx = x * t1 / r;
                dt2_dy = y * t1 / r;
                dt2_dz = z * t1 / r;

                beptr1x[idx].re = cc * (t2 * dy1_dx + y1 * dt2_dx);
                beptr1y[idx].re = cc * (t2 * dy1_dy + y1 * dt2_dy);
                beptr1z[idx].re = cc * (t2 * dy1_dz + y1 * dt2_dz);
                beptr1x[idx].im = 0.0;
                beptr1y[idx].im = 0.0;
                beptr1z[idx].im = 0.0;

                beptr2x[idx].re = cc * (t2 * dy2_dx + y2 * dt2_dx);
                beptr2y[idx].re = cc * (t2 * dy2_dy + y2 * dt2_dy);
                beptr2z[idx].re = cc * (t2 * dy2_dz + y2 * dt2_dz);
                beptr2x[idx].im = 0.0;
                beptr2y[idx].im = 0.0;
                beptr2z[idx].im = 0.0;

                beptr3x[idx].re = cc * (t2 * dy3_dx + y3 * dt2_dx);
                beptr3y[idx].re = cc * (t2 * dy3_dy + y3 * dt2_dy);
                beptr3z[idx].re = cc * (t2 * dy3_dz + y3 * dt2_dz);
                beptr3x[idx].im = 0.0;
                beptr3y[idx].im = 0.0;
                beptr3z[idx].im = 0.0;

                beptr4x[idx].re = cc * (t2 * dy4_dx + y4 * dt2_dx);
                beptr4y[idx].re = cc * (t2 * dy4_dy + y4 * dt2_dy);
                beptr4z[idx].re = cc * (t2 * dy4_dz + y4 * dt2_dz);
                beptr4x[idx].im = 0.0;
                beptr4y[idx].im = 0.0;
                beptr4z[idx].im = 0.0;

                beptr5x[idx].re = cc * (t2 * dy5_dx + y5 * dt2_dx);
                beptr5y[idx].re = cc * (t2 * dy5_dy + y5 * dt2_dy);
                beptr5z[idx].re = cc * (t2 * dy5_dz + y5 * dt2_dz);
                beptr5x[idx].im = 0.0;
                beptr5y[idx].im = 0.0;
                beptr5z[idx].im = 0.0;

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

    fftwnd_one(p1, beptr1y, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betay1);

    fftwnd_one(p1, beptr1z, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betaz1);

    fftwnd_one(p1, beptr2x, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betax2);

    fftwnd_one(p1, beptr2y, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betay2);

    fftwnd_one(p1, beptr2z, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betaz2);

    fftwnd_one(p1, beptr3x, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betax3);

    fftwnd_one(p1, beptr3y, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betay3);

    fftwnd_one(p1, beptr3z, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betaz3);

    fftwnd_one(p1, beptr4x, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betax4);

    fftwnd_one(p1, beptr4y, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betay4);

    fftwnd_one(p1, beptr4z, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betaz4);

    fftwnd_one(p1, beptr5x, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betax5);

    fftwnd_one(p1, beptr5y, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betay5);

    fftwnd_one(p1, beptr5z, gbptr);
    pack_gftoc(sp, gbptr, gwptr);
    fftwnd_one(p2, gwptr, weptr);
    assign_weight(sp, weptr, betaz5);


    my_free(beptr1x);
    my_free(beptr2x);
    my_free(beptr3x);
    my_free(beptr4x);
    my_free(beptr5x);
    my_free(gbptr);
    my_free(weptr);
    my_free(gwptr);

}
