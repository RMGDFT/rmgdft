/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"

void init_derweight_d (SPECIES * sp,
                       fftw_complex * rtptr_x,
                       fftw_complex * rtptr_y, fftw_complex * rtptr_z, int ip, fftwnd_plan p1)
{
#if !FDIFF_BETA

    int idx, ix, iy, iz, size, coarse_size, iend, ibegin;
    REAL r, ax[3], bx[3], xc, yc, zc, t1, t2, t3, invdr;
    REAL x, y, z, rsq, rsqd, r4, xsq, ysq, zsq, y1, y2, y3, y4, y5;
    REAL dy1_dx, dy1_dy, dy1_dz, dy2_dx, dy2_dy, dy2_dz, dy3_dx, dy3_dy, dy3_dz, dy4_dx, dy4_dy,
        dy4_dz, dy5_dx, dy5_dy, dy5_dz;
    REAL dt2_dx, dt2_dy, dt2_dz;
    REAL cc, hxx, hyy, hzz;
    fftw_complex *weptr1x, *weptr1y, *weptr1z, *gwptr;
    fftw_complex *weptr2x, *weptr2y, *weptr2z;
    fftw_complex *weptr3x, *weptr3y, *weptr3z;
    fftw_complex *weptr4x, *weptr4y, *weptr4z;
    fftw_complex *weptr5x, *weptr5y, *weptr5z;
    fftw_complex *r1x, *r1y, *r1z;
    fftw_complex *r2x, *r2y, *r2z;
    fftw_complex *r3x, *r3y, *r3z;
    fftw_complex *r4x, *r4y, *r4z;
    fftw_complex *r5x, *r5y, *r5z;

    invdr = 1.0 / sp->drnlig;


    /*Number of grid points in the non-local box in coarse and double grids */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    my_malloc (weptr1x, 16 * size, fftw_complex);
    if (weptr1x == NULL)
        error_handler ("can't allocate memory\n");

    weptr1y = weptr1x + size;
    weptr1z = weptr1y + size;

    weptr2x = weptr1z + size;
    weptr2y = weptr2x + size;
    weptr2z = weptr2y + size;

    weptr3x = weptr2z + size;
    weptr3y = weptr3x + size;
    weptr3z = weptr3y + size;

    weptr4x = weptr3z + size;
    weptr4y = weptr4x + size;
    weptr4z = weptr4y + size;

    weptr5x = weptr4z + size;
    weptr5y = weptr5x + size;
    weptr5z = weptr5y + size;


    gwptr = weptr5z + size;


    hxx = ct.hxgrid / (REAL) ct.nxfgrid;
    hyy = ct.hygrid / (REAL) ct.nyfgrid;
    hzz = ct.hzgrid / (REAL) ct.nzfgrid;

    r1x = rtptr_x;
    r1y = rtptr_y;
    r1z = rtptr_z;

    r2x = r1x + coarse_size;
    r2y = r1y + coarse_size;
    r2z = r1z + coarse_size;

    r3x = r2x + coarse_size;
    r3y = r2y + coarse_size;
    r3z = r2z + coarse_size;

    r4x = r3x + coarse_size;
    r4y = r3y + coarse_size;
    r4z = r3z + coarse_size;

    r5x = r4x + coarse_size;
    r5y = r4y + coarse_size;
    r5z = r4z + coarse_size;


    cc = sqrt (5.0 / (4.0 * PI));
    t3 = sqrt (3.0);

    ibegin = -(sp->nldim / 2) * ct.nxfgrid;
    iend = ibegin + sp->nlfdim;

    idx = 0;
    for (ix = ibegin; ix < iend; ix++)
    {
        xc = (REAL) ix *hxx;

        for (iy = ibegin; iy < iend; iy++)
        {
            yc = (REAL) iy *hyy;

            for (iz = ibegin; iz < iend; iz++)
            {
                zc = (REAL) iz *hzz;


                ax[0] = xc;
                ax[1] = yc;
                ax[2] = zc;

                to_cartesian (ax, bx);
                r = metric (ax);
                x = bx[0];
                y = bx[1];
                z = bx[2];

                rsq = x * x + y * y + z * z;
                /*This is apparently not used */
                /*r3 = r * rsq + 1.0e-30; */
                rsqd = rsq + 1.0e-20;
                r4 = rsq * rsq + 1.e-40;

                t1 = t3 * linint (&sp->drbetalig[ip][0], r, invdr);
                t2 = t3 * linint (&sp->betalig[ip][0], r, invdr);
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


                weptr1x[idx].re = cc * (t2 * dy1_dx + y1 * dt2_dx);
                weptr1y[idx].re = cc * (t2 * dy1_dy + y1 * dt2_dy);
                weptr1z[idx].re = cc * (t2 * dy1_dz + y1 * dt2_dz);
                weptr1x[idx].im = 0.0;
                weptr1y[idx].im = 0.0;
                weptr1z[idx].im = 0.0;

                weptr2x[idx].re = cc * (t2 * dy2_dx + y2 * dt2_dx);
                weptr2y[idx].re = cc * (t2 * dy2_dy + y2 * dt2_dy);
                weptr2z[idx].re = cc * (t2 * dy2_dz + y2 * dt2_dz);
                weptr2x[idx].im = 0.0;
                weptr2y[idx].im = 0.0;
                weptr2z[idx].im = 0.0;

                weptr3x[idx].re = cc * (t2 * dy3_dx + y3 * dt2_dx);
                weptr3y[idx].re = cc * (t2 * dy3_dy + y3 * dt2_dy);
                weptr3z[idx].re = cc * (t2 * dy3_dz + y3 * dt2_dz);
                weptr3x[idx].im = 0.0;
                weptr3y[idx].im = 0.0;
                weptr3z[idx].im = 0.0;

                weptr4x[idx].re = cc * (t2 * dy4_dx + y4 * dt2_dx);
                weptr4y[idx].re = cc * (t2 * dy4_dy + y4 * dt2_dy);
                weptr4z[idx].re = cc * (t2 * dy4_dz + y4 * dt2_dz);
                weptr4x[idx].im = 0.0;
                weptr4y[idx].im = 0.0;
                weptr4z[idx].im = 0.0;

                weptr5x[idx].re = cc * (t2 * dy5_dx + y5 * dt2_dx);
                weptr5y[idx].re = cc * (t2 * dy5_dy + y5 * dt2_dy);
                weptr5z[idx].re = cc * (t2 * dy5_dz + y5 * dt2_dz);
                weptr5x[idx].im = 0.0;
                weptr5y[idx].im = 0.0;
                weptr5z[idx].im = 0.0;


                idx++;
            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    fftwnd_one (p1, weptr1x, gwptr);
    pack_gftoc (sp, gwptr, r1x);

    fftwnd_one (p1, weptr1y, gwptr);
    pack_gftoc (sp, gwptr, r1y);

    fftwnd_one (p1, weptr1z, gwptr);
    pack_gftoc (sp, gwptr, r1z);



    fftwnd_one (p1, weptr2x, gwptr);
    pack_gftoc (sp, gwptr, r2x);

    fftwnd_one (p1, weptr2y, gwptr);
    pack_gftoc (sp, gwptr, r2y);

    fftwnd_one (p1, weptr2z, gwptr);
    pack_gftoc (sp, gwptr, r2z);



    fftwnd_one (p1, weptr3x, gwptr);
    pack_gftoc (sp, gwptr, r3x);

    fftwnd_one (p1, weptr3y, gwptr);
    pack_gftoc (sp, gwptr, r3y);

    fftwnd_one (p1, weptr3z, gwptr);
    pack_gftoc (sp, gwptr, r3z);



    fftwnd_one (p1, weptr4x, gwptr);
    pack_gftoc (sp, gwptr, r4x);

    fftwnd_one (p1, weptr4y, gwptr);
    pack_gftoc (sp, gwptr, r4y);

    fftwnd_one (p1, weptr4z, gwptr);
    pack_gftoc (sp, gwptr, r4z);



    fftwnd_one (p1, weptr5x, gwptr);
    pack_gftoc (sp, gwptr, r5x);

    fftwnd_one (p1, weptr5y, gwptr);
    pack_gftoc (sp, gwptr, r5y);

    fftwnd_one (p1, weptr5z, gwptr);
    pack_gftoc (sp, gwptr, r5z);





    my_free (weptr1x);
#endif
}
