/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "grid.h"
#include "const.h"
#include "params.h"
#include "rmgtypes.h"
#include "rmg_alloc.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"

void init_derweight_d (SPECIES * sp,
                       double complex * rtptr_x,
                       double complex * rtptr_y, fftw_complex * rtptr_z, int ip, fftw_plan p1)
{
#if !FDIFF_BETA

    int idx, ix, iy, iz, size, coarse_size, iend, ibegin;
    rmg_double_t r, ax[3], bx[3], xc, yc, zc, t1, t2, t3, invdr;
    rmg_double_t x, y, z, rsq, rsqd, r4, xsq, ysq, zsq, y1, y2, y3, y4, y5;
    rmg_double_t dy1_dx, dy1_dy, dy1_dz, dy2_dx, dy2_dy, dy2_dz, dy3_dx, dy3_dy, dy3_dz, dy4_dx, dy4_dy,
        dy4_dz, dy5_dx, dy5_dy, dy5_dz;
    rmg_double_t dt2_dx, dt2_dy, dt2_dz;
    rmg_double_t cc, hxx, hyy, hzz;
    double complex *weptr1x, *weptr1y, *weptr1z, *gwptr;
    double complex *weptr2x, *weptr2y, *weptr2z;
    double complex *weptr3x, *weptr3y, *weptr3z;
    double complex *weptr4x, *weptr4y, *weptr4z;
    double complex *weptr5x, *weptr5y, *weptr5z;
    double complex *r1x, *r1y, *r1z;
    double complex *r2x, *r2y, *r2z;
    double complex *r3x, *r3y, *r3z;
    double complex *r4x, *r4y, *r4z;
    double complex *r5x, *r5y, *r5z;

    invdr = 1.0 / sp->drnlig;


    /*Number of grid points in the non-local box in coarse and double grids */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    my_malloc (weptr1x, 16 * size, double complex);
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


    hxx = get_hxgrid() / (rmg_double_t) ct.nxfgrid;
    hyy = get_hygrid() / (rmg_double_t) ct.nyfgrid;
    hzz = get_hzgrid() / (rmg_double_t) ct.nzfgrid;

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
        xc = (rmg_double_t) ix *hxx;

        for (iy = ibegin; iy < iend; iy++)
        {
            yc = (rmg_double_t) iy *hyy;

            for (iz = ibegin; iz < iend; iz++)
            {
                zc = (rmg_double_t) iz *hzz;


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


                weptr1x[idx] = cc * (t2 * dy1_dx + y1 * dt2_dx) + 0.0I;
                weptr1y[idx] = cc * (t2 * dy1_dy + y1 * dt2_dy) + 0.0I;
                weptr1z[idx] = cc * (t2 * dy1_dz + y1 * dt2_dz) + 0.0I;

                weptr2x[idx] = cc * (t2 * dy2_dx + y2 * dt2_dx) + 0.0I;
                weptr2y[idx] = cc * (t2 * dy2_dy + y2 * dt2_dy) + 0.0I;
                weptr2z[idx] = cc * (t2 * dy2_dz + y2 * dt2_dz) + 0.0I;

                weptr3x[idx] = cc * (t2 * dy3_dx + y3 * dt2_dx) + 0.0I;
                weptr3y[idx] = cc * (t2 * dy3_dy + y3 * dt2_dy) + 0.0I;
                weptr3z[idx] = cc * (t2 * dy3_dz + y3 * dt2_dz) + 0.0I;

                weptr4x[idx] = cc * (t2 * dy4_dx + y4 * dt2_dx) + 0.0I;
                weptr4y[idx] = cc * (t2 * dy4_dy + y4 * dt2_dy) + 0.0I;
                weptr4z[idx] = cc * (t2 * dy4_dz + y4 * dt2_dz) + 0.0I;

                weptr5x[idx] = cc * (t2 * dy5_dx + y5 * dt2_dx) + 0.0I;
                weptr5y[idx] = cc * (t2 * dy5_dy + y5 * dt2_dy) + 0.0I;
                weptr5z[idx] = cc * (t2 * dy5_dz + y5 * dt2_dz) + 0.0I;


                idx++;
            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    fftw_execute_dft (p1, weptr1x, gwptr);
    pack_gftoc (sp, gwptr, r1x);

    fftw_execute_dft (p1, weptr1y, gwptr);
    pack_gftoc (sp, gwptr, r1y);

    fftw_execute_dft (p1, weptr1z, gwptr);
    pack_gftoc (sp, gwptr, r1z);



    fftw_execute_dft (p1, weptr2x, gwptr);
    pack_gftoc (sp, gwptr, r2x);

    fftw_execute_dft (p1, weptr2y, gwptr);
    pack_gftoc (sp, gwptr, r2y);

    fftw_execute_dft (p1, weptr2z, gwptr);
    pack_gftoc (sp, gwptr, r2z);



    fftw_execute_dft (p1, weptr3x, gwptr);
    pack_gftoc (sp, gwptr, r3x);

    fftw_execute_dft (p1, weptr3y, gwptr);
    pack_gftoc (sp, gwptr, r3y);

    fftw_execute_dft (p1, weptr3z, gwptr);
    pack_gftoc (sp, gwptr, r3z);



    fftw_execute_dft (p1, weptr4x, gwptr);
    pack_gftoc (sp, gwptr, r4x);

    fftw_execute_dft (p1, weptr4y, gwptr);
    pack_gftoc (sp, gwptr, r4y);

    fftw_execute_dft (p1, weptr4z, gwptr);
    pack_gftoc (sp, gwptr, r4z);



    fftw_execute_dft (p1, weptr5x, gwptr);
    pack_gftoc (sp, gwptr, r5x);

    fftw_execute_dft (p1, weptr5y, gwptr);
    pack_gftoc (sp, gwptr, r5y);

    fftw_execute_dft (p1, weptr5z, gwptr);
    pack_gftoc (sp, gwptr, r5z);





    my_free (weptr1x);
#endif
}
