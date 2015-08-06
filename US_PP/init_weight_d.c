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

void init_weight_d (SPECIES * sp, fftw_complex * rtptr, int ip, fftw_plan p1)
{

    int idx, ix, iy, iz, size, coarse_size, iend, ibegin;
    double r, ax[3], bx[3], xc, yc, zc, t1, t2, rsq1, invdr;
    double cc, hxx, hyy, hzz;
    double complex *weptr1, *weptr2, *weptr3, *weptr4, *weptr5, *gwptr;
    double complex *r1, *r2, *r3, *r4, *r5;
    int ixx, iyy, izz;

    invdr = 1.0 / sp->drnlig;


    /*Number of grid points in th enon-local box in coarse and double grids */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    weptr1 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr2 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr3 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr4 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr5 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    gwptr = (double complex *)fftw_malloc(sizeof(double complex) * size);

    if ((weptr1 == NULL) || (weptr2 == NULL) || (weptr3 == NULL) || (weptr4 == NULL) || (weptr5 == NULL) || (gwptr == NULL))
        error_handler ("can't allocate memory\n");

    hxx = get_hxgrid() / (double) ct.nxfgrid;
    hyy = get_hygrid() / (double) ct.nyfgrid;
    hzz = get_hzgrid() / (double) ct.nzfgrid;

    r1 = rtptr;
    r2 = r1 + coarse_size;
    r3 = r2 + coarse_size;
    r4 = r3 + coarse_size;
    r5 = r4 + coarse_size;

    cc = sqrt (5.0 / (4.0 * PI));
    t2 = sqrt (3.0);

    ibegin = -sp->nlfdim / 2;
    iend = ibegin + sp->nlfdim;

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

                to_cartesian (ax, bx);
                r = metric (ax);
                rsq1 = r * r + 1.0e-20;
                t1 = linint (&sp->betalig[ip][0], r, invdr);
                t1 = t1 * t2;
                weptr1[idx] = cc * t1 * bx[0] * bx[1] / rsq1 + 0.0I;
                weptr2[idx] = cc * t1 * bx[0] * bx[2] / rsq1 + 0.0I;
                weptr3[idx] = cc * t1 * (t2 * bx[2] * bx[2] - rsq1 / t2) / (2.0 * rsq1) + 0.0I;
                weptr4[idx] = cc * t1 * bx[1] * bx[2] / rsq1 + 0.0I;
                weptr5[idx] = cc * t1 * (bx[0] * bx[0] - bx[1] * bx[1]) / (2.0 * rsq1) + 0.0I;


                if((ix*2 + sp->nlfdim) == 0 | (iy*2 + sp->nlfdim) == 0 | (iz*2 + sp->nlfdim) == 0 )
    {
                    weptr1[idx] = 0.0;
                    weptr2[idx] = 0.0;
                    weptr3[idx] = 0.0;
                    weptr4[idx] = 0.0;
                    weptr5[idx] = 0.0;
}
            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    fftw_execute_dft (p1, weptr1, gwptr);
    pack_gftoc (sp, gwptr, r1);

    fftw_execute_dft (p1, weptr2, gwptr);
    pack_gftoc (sp, gwptr, r2);

    fftw_execute_dft (p1, weptr3, gwptr);
    pack_gftoc (sp, gwptr, r3);

    fftw_execute_dft (p1, weptr4, gwptr);
    pack_gftoc (sp, gwptr, r4);

    fftw_execute_dft (p1, weptr5, gwptr);
    pack_gftoc (sp, gwptr, r5);

    fftw_free (gwptr);
    fftw_free (weptr5);
    fftw_free (weptr4);
    fftw_free (weptr3);
    fftw_free (weptr2);
    fftw_free (weptr1);

}
