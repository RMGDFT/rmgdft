/************************** SVN Revision Information **************************
 **    $Id: init_weight_p.c 1066 2009-08-31 18:41:09Z froze $    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"

void init_weight_p (SPECIES * sp, fftw_complex * rtptr, int ip, fftwnd_plan p1)
{

    int idx, ix, iy, iz, size, coarse_size, ibegin, iend;
    REAL r, ax[3], bx[3], xc, yc, zc, cc, t1, invdr;
    REAL hxx, hyy, hzz;
    fftw_complex *weptr1, *weptr2, *weptr3, *gwptr;
    fftw_complex *r1, *r2, *r3;

    invdr = 1.0 / sp->drnlig;


    /*Number of grid points in th enon-local box in coarse and double grids */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    my_malloc (weptr1, 4 * size, fftw_complex);
    if (weptr1 == NULL)
        error_handler ("can't allocate memory\n");

    weptr2 = weptr1 + size;
    weptr3 = weptr2 + size;
    gwptr = weptr3 + size;

    hxx = ct.hxgrid / (REAL) ct.nxfgrid;
    hyy = ct.hygrid / (REAL) ct.nyfgrid;
    hzz = ct.hzgrid / (REAL) ct.nzfgrid;

    r1 = rtptr;
    r2 = r1 + coarse_size;
    r3 = r2 + coarse_size;

    cc = sqrt (3.0 / (4.0 * PI));

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
                t1 = linint (&sp->betalig[ip][0], r, invdr);
                to_cartesian (ax, bx);
                r += 1.0e-10;

                weptr1[idx].re = cc * bx[0] * t1 / r;
                weptr2[idx].re = cc * bx[2] * t1 / r;
                weptr3[idx].re = cc * bx[1] * t1 / r;
                weptr1[idx].im = 0.0;
                weptr2[idx].im = 0.0;
                weptr3[idx].im = 0.0;

                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    fftwnd_one (p1, weptr1, gwptr);
    pack_gftoc (sp, gwptr, r1);

    fftwnd_one (p1, weptr2, gwptr);
    pack_gftoc (sp, gwptr, r2);

    fftwnd_one (p1, weptr3, gwptr);
    pack_gftoc (sp, gwptr, r3);

    my_free (weptr1);

}
