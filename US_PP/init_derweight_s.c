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
#include "rmg_alloc.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"

void init_derweight_s (SPECIES * sp,
                       fftw_complex * rtptr_x,
                       fftw_complex * rtptr_y, fftw_complex * rtptr_z, int ip, fftw_plan p1)
{

    int idx, ix, iy, iz, size, ibegin, iend;
    double r, ax[3], bx[3], xc, yc, zc;
    double invdr, t1, hxx, hyy, hzz;
    double complex *weptrx, *weptry, *weptrz, *gwptr;
    int ixx, iyy, izz;



    /* nlffdim is size of the non-local box in the double grid */
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    my_malloc (weptrx, 4 * size, fftw_complex);
    if (weptrx == NULL)
        error_handler ("can't allocate memory\n");
    weptry = weptrx + size;
    weptrz = weptry + size;
    gwptr = weptrz + size;

    hxx = get_hxgrid() / (double) ct.nxfgrid;
    hyy = get_hygrid() / (double) ct.nyfgrid;
    hzz = get_hzgrid() / (double) ct.nzfgrid;

    invdr = 1.0 / sp->drnlig;

    /*We assume that ion is in the center of non-local box */
    ibegin = -sp->nlfdim / 2;
    iend = ibegin + sp->nlfdim;

    idx = 0;
    for (ix = ibegin; ix < iend; ix++)
    {
        ixx = ix;
        if (ixx < 0) ixx = ix + sp->nlfdim;
        xc = (double) ix *hxx;

        for (iy = ibegin; iy < iend; iy++)
        {
            iyy = iy;
            if (iyy < 0) iyy = iy + sp->nlfdim;
            yc = (double) iy *hyy;

            for (iz = ibegin; iz < iend; iz++)
            {

                izz = iz;
                if (izz < 0) izz = iz + sp->nlfdim;
                idx = ixx *sp->nlfdim * sp->nlfdim + iyy * sp->nlfdim + izz;

                zc = (double) iz *hzz;


                ax[0] = xc;
                ax[1] = yc;
                ax[2] = zc;

                r = metric (ax);
                to_cartesian (ax, bx);
                r += 1.0e-10;

                t1 = linint (&sp->drbetalig[ip][0], r, invdr);


                weptrx[idx] = sqrt (1.0 / (4.0 * PI)) * t1 * bx[0] / r + 0.0I;
                weptry[idx] = sqrt (1.0 / (4.0 * PI)) * t1 * bx[1] / r + 0.0I;
                weptrz[idx] = sqrt (1.0 / (4.0 * PI)) * t1 * bx[2] / r + 0.0I;


                if((ix*2 + sp->nlfdim) == 0 | (iy*2 + sp->nlfdim) == 0 | (iz*2 + sp->nlfdim) == 0 )
                {
                    weptrx[idx] = 0.0;
                    weptry[idx] = 0.0;
                    weptrz[idx] = 0.0;
                }

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    /*Fourier transform and restricting G-space for all three derivatives */
    fftw_execute_dft (p1, weptrx, gwptr);
    pack_gftoc (sp, gwptr, rtptr_x);

    fftw_execute_dft (p1, weptry, gwptr);
    pack_gftoc (sp, gwptr, rtptr_y);

    fftw_execute_dft (p1, weptrz, gwptr);
    pack_gftoc (sp, gwptr, rtptr_z);

    my_free (weptrx);
}
