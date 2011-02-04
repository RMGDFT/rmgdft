/************************** SVN Revision Information **************************
 **    $Id: init_weight_d.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"

void init_weight_d (SPECIES *sp, fftw_complex *rtptr, int ip, fftwnd_plan p1)
{

    int idx, ix, iy, iz, size, coarse_size, iend, ibegin;
    REAL r, ax[3], bx[3], xc, yc, zc, t1, t2, rsq1, invdr;
    REAL cc, hxx, hyy, hzz;
    fftw_complex *weptr1, *weptr2, *weptr3, *weptr4, *weptr5, *gwptr;
    fftw_complex *r1, *r2, *r3, *r4, *r5;

    invdr = 1.0 / sp->drnlig;


    /*Number of grid points in th enon-local box in coarse and double grids */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    my_malloc( weptr1, 6 * size, fftw_complex );
    if (weptr1 == NULL)
        error_handler ("can't allocate memory\n");

    weptr2 = weptr1 + size;
    weptr3 = weptr2 + size;
    weptr4 = weptr3 + size;
    weptr5 = weptr4 + size;
    gwptr = weptr5 + size;

    hxx = ct.hxgrid / (REAL) ct.nxfgrid;
    hyy = ct.hygrid / (REAL) ct.nyfgrid;
    hzz = ct.hzgrid / (REAL) ct.nzfgrid;

    r1 = rtptr;
    r2 = r1 + coarse_size;
    r3 = r2 + coarse_size;
    r4 = r3 + coarse_size;
    r5 = r4 + coarse_size;

    cc = sqrt (5.0 / (4.0 * PI));
    t2 = sqrt (3.0);

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
                rsq1 = r * r + 1.0e-20;
                t1 = linint (&sp->betalig[ip][0], r, invdr);
                t1 = t1 * t2;
                weptr1[idx].re = cc * t1 * bx[0] * bx[1] / rsq1;
                weptr2[idx].re = cc * t1 * bx[0] * bx[2] / rsq1;
                weptr3[idx].re = cc * t1 * (t2 * bx[2] * bx[2] - rsq1 / t2) / (2.0 * rsq1);
                weptr4[idx].re = cc * t1 * bx[1] * bx[2] / rsq1;
                weptr5[idx].re = cc * t1 * (bx[0] * bx[0] - bx[1] * bx[1]) / (2.0 * rsq1);
                weptr1[idx].im = 0.0;
                weptr2[idx].im = 0.0;
                weptr3[idx].im = 0.0;
                weptr4[idx].im = 0.0;
                weptr5[idx].im = 0.0;

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

    fftwnd_one (p1, weptr4, gwptr);
    pack_gftoc (sp, gwptr, r4);

    fftwnd_one (p1, weptr5, gwptr);
    pack_gftoc (sp, gwptr, r5);

    my_free(weptr1);

}
