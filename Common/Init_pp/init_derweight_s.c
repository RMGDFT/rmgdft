/************************** SVN Revision Information **************************
 **    $Id: init_derweight_s.c 1066 2009-08-31 18:41:09Z froze $    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"

void init_derweight_s (SPECIES * sp,
                       fftw_complex * rtptr_x,
                       fftw_complex * rtptr_y, fftw_complex * rtptr_z, int ip, fftwnd_plan p1)
{
#if !FDIFF_BETA

    int idx, ix, iy, iz, size, ibegin, iend;
    REAL r, ax[3], bx[3], xc, yc, zc;
    REAL invdr, t1, hxx, hyy, hzz;
    fftw_complex *weptrx, *weptry, *weptrz, *gwptr;



    /* nlffdim is size of the non-local box in the double grid */
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    my_malloc (weptrx, 4 * size, fftw_complex);
    if (weptrx == NULL)
        error_handler ("can't allocate memory\n");
    weptry = weptrx + size;
    weptrz = weptry + size;
    gwptr = weptrz + size;

    hxx = ct.hxgrid / (REAL) ct.nxfgrid;
    hyy = ct.hygrid / (REAL) ct.nyfgrid;
    hzz = ct.hzgrid / (REAL) ct.nzfgrid;

    invdr = 1.0 / sp->drnlig;

    /*We assume that ion is in the center of non-local box */
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

                r = metric (ax);
                to_cartesian (ax, bx);
                r += 1.0e-10;

                t1 = linint (&sp->drbetalig[ip][0], r, invdr);


                weptrx[idx].re = sqrt (1.0 / (4.0 * PI)) * t1 * bx[0] / r;
                weptry[idx].re = sqrt (1.0 / (4.0 * PI)) * t1 * bx[1] / r;
                weptrz[idx].re = sqrt (1.0 / (4.0 * PI)) * t1 * bx[2] / r;
                weptrx[idx].im = 0.0;
                weptry[idx].im = 0.0;
                weptrz[idx].im = 0.0;

                idx++;
            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    /*Fourier transform and restricting G-space for all three derivatives */
    fftwnd_one (p1, weptrx, gwptr);
    pack_gftoc (sp, gwptr, rtptr_x);

    fftwnd_one (p1, weptry, gwptr);
    pack_gftoc (sp, gwptr, rtptr_y);

    fftwnd_one (p1, weptrz, gwptr);
    pack_gftoc (sp, gwptr, rtptr_z);

    my_free (weptrx);
#endif
}
