/************************** SVN Revision Information **************************
 **    $Id: init_weight_s.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"

void init_weight_s (SPECIES *sp, fftw_complex *rtptr, int ip, fftwnd_plan p1)
{

    int idx, ix, iy, iz, size, ibegin, iend;
    REAL r, ax[3], xc, yc, zc;
    REAL invdr, t1, hxx, hyy, hzz;
    fftw_complex *weptr, *gwptr;



    /* nlfdim is size of the non-local box in the double grid */
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    my_malloc( weptr, 2 * size, fftw_complex );
    if (weptr == NULL)
        error_handler ("can't allocate memory\n");
    gwptr = weptr + size;

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
                t1 = linint (&sp->betalig[ip][0], r, invdr);
                weptr[idx].re = sqrt (1.0 / (4.0 * PI)) * t1;
                weptr[idx].im = 0.0;

                idx++;
            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    fftwnd_one (p1, weptr, gwptr);
    pack_gftoc (sp, gwptr, rtptr);

    my_free(weptr);
}
