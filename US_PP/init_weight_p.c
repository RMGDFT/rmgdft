/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "main.h"
#include "common_prototypes.h"

void init_weight_p (SPECIES * sp, fftw_complex * rtptr, int ip, fftw_plan p1)
{

    int idx, ix, iy, iz, size, coarse_size, ibegin, iend;
    rmg_double_t r, ax[3], bx[3], xc, yc, zc, cc, t1, invdr;
    rmg_double_t hxx, hyy, hzz;
    double complex *weptr1, *weptr2, *weptr3, *gwptr;
    double complex *r1, *r2, *r3;
    int ixx, iyy, izz;
    double complex *in, *out;

    invdr = 1.0 / sp->drnlig;


        /*This is something we need to do only once per species, so do not use wisdom */
//    in = (double complex *)fftw_malloc(sizeof(double complex) * sp->nlfdim * sp->nlfdim * sp->nlfdim);
//    out = (double complex *)fftw_malloc(sizeof(double complex) * sp->nlfdim * sp->nlfdim * sp->nlfdim);


 //   p1 = fftw_plan_dft_3d (sp->nlfdim, sp->nlfdim, sp->nlfdim, in, out, FFTW_FORWARD, FFTW_MEASURE);


    /*Number of grid points in th enon-local box in coarse and double grids */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    weptr1 = (double complex *)fftw_malloc(sizeof(double complex) * 4 * size);

    if (weptr1 == NULL)
        error_handler ("can't allocate memory\n");

    weptr2 = weptr1 + size;
    weptr3 = weptr2 + size;
    gwptr = weptr3 + size;

    hxx = get_hxgrid() / (rmg_double_t) ct.nxfgrid;
    hyy = get_hygrid() / (rmg_double_t) ct.nyfgrid;
    hzz = get_hzgrid() / (rmg_double_t) ct.nzfgrid;

    r1 = rtptr;
    r2 = r1 + coarse_size;
    r3 = r2 + coarse_size;

    cc = sqrt (3.0 / (4.0 * PI));

    ibegin = -sp->nlfdim/2;
    iend = ibegin + sp->nlfdim;

    for (ix = ibegin; ix < iend; ix++)
    {
        ixx = ix;
        if (ixx < 0) ixx = ix + sp->nlfdim;

        xc = (rmg_double_t) ix *hxx;

        for (iy = ibegin; iy < iend; iy++)
        {
            iyy = iy;
            if (iyy < 0) iyy = iy + sp->nlfdim;

            yc = (rmg_double_t) iy *hyy;

            for (iz = ibegin; iz < iend; iz++)
            {
                izz = iz;
                if (izz < 0) izz = iz + sp->nlfdim;
                idx = ixx *sp->nlfdim * sp->nlfdim + iyy * sp->nlfdim + izz;

                zc = (rmg_double_t) iz *hzz;

                ax[0] = xc;
                ax[1] = yc;
                ax[2] = zc;

                r = metric (ax);
                t1 = linint (&sp->betalig[ip][0], r, invdr);
                to_cartesian (ax, bx);
                r += 1.0e-10;

                weptr1[idx] = cc * bx[0] * t1 / r + 0.0I;
                weptr2[idx] = cc * bx[2] * t1 / r + 0.0I;
                weptr3[idx] = cc * bx[1] * t1 / r + 0.0I;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    fftw_execute_dft (p1, weptr1, gwptr);
    pack_gftoc (sp, gwptr, r1);

    fftw_execute_dft (p1, weptr2, gwptr);
    pack_gftoc (sp, gwptr, r2);


    fftw_execute_dft (p1, weptr3, gwptr);
    pack_gftoc (sp, gwptr, r3);

    fftw_free (weptr1);
    fftw_free (in);
    fftw_free (out);

}
