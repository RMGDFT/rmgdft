/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "const.h"
#include "params.h"
#include "rmgtypes.h"
#include "rmg_alloc.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "prototypes.h"

void init_derweight_p (SPECIES * sp,
                       fftw_complex * rtptr_x,
                       fftw_complex * rtptr_y, fftw_complex * rtptr_z, int ip, fftwnd_plan p1)
{
#if !FDIFF_BETA

    int idx, ix, iy, iz, size, coarse_size, ibegin, iend;
    rmg_double_t r, rsq, r3, rsqd, ax[3], bx[3], x, y, z, xc, yc, zc, cc, t1, t2, invdr;
    rmg_double_t hxx, hyy, hzz;
    fftw_complex *weptr1x, *weptr1y, *weptr1z, *gwptr;
    fftw_complex *weptr2x, *weptr2y, *weptr2z;
    fftw_complex *weptr3x, *weptr3y, *weptr3z;
    fftw_complex *r1x, *r1y, *r1z, *r2x, *r2y, *r2z, *r3x, *r3y, *r3z;

    invdr = 1.0 / sp->drnlig;


    /*Number of grid points in th enon-local box in coarse and double grids */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    my_malloc (weptr1x, 10 * size, fftw_complex);
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

    gwptr = weptr3z + size;

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

    cc = sqrt (3.0 / (4.0 * PI));

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

                r = metric (ax);
                t1 = linint (&sp->drbetalig[ip][0], r, invdr);
                t2 = linint (&sp->betalig[ip][0], r, invdr);

                to_cartesian (ax, bx);
                x = bx[0];
                y = bx[1];
                z = bx[2];

                rsq = x * x + y * y + z * z;
                r3 = r * rsq + 1.0e-20;
                rsqd = rsq + 1.0e-20;


                weptr1x[idx].re = cc * ((t2 * (rsq - x * x) / r3) + (x * x * t1 / rsqd));
                weptr1y[idx].re = cc * ((-t2 * x * y / r3) + (x * y * t1 / rsqd));
                weptr1z[idx].re = cc * ((-t2 * x * z / r3) + (x * z * t1 / rsqd));;

                weptr2x[idx].re = cc * ((-t2 * x * z / r3) + (x * z * t1 / rsqd));
                weptr2y[idx].re = cc * ((-t2 * y * z / r3) + (y * z * t1 / rsqd));
                weptr2z[idx].re = cc * ((t2 * (rsq - z * z) / r3) + (z * z * t1 / rsqd));

                weptr3x[idx].re = cc * ((-t2 * x * y / r3) + (x * y * t1 / rsqd));
                weptr3y[idx].re = cc * ((t2 * (rsq - y * y) / r3) + (y * y * t1 / rsqd));
                weptr3z[idx].re = cc * ((-t2 * y * z / r3) + (y * z * t1 / rsqd));




                weptr1x[idx].im = 0.0;
                weptr1y[idx].im = 0.0;
                weptr1z[idx].im = 0.0;

                weptr2x[idx].im = 0.0;
                weptr2y[idx].im = 0.0;
                weptr2z[idx].im = 0.0;

                weptr3x[idx].im = 0.0;
                weptr3y[idx].im = 0.0;
                weptr3z[idx].im = 0.0;

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



    my_free (weptr1x);

#endif
}
