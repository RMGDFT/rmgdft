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
    double r, ax[3], bx[3], xc, yc, zc, cc, t1;
    double hxx, hyy, hzz;
    double complex *weptr1, *weptr2, *weptr3, *gwptr;
    double complex *r1, *r2, *r3;
    int ixx, iyy, izz;


        /*This is something we need to do only once per species, so do not use wisdom */
//    in = (double complex *)fftw_malloc(sizeof(double complex) * sp->nlfdim * sp->nlfdim * sp->nlfdim);
//    out = (double complex *)fftw_malloc(sizeof(double complex) * sp->nlfdim * sp->nlfdim * sp->nlfdim);


 //   p1 = fftw_plan_dft_3d (sp->nlfdim, sp->nlfdim, sp->nlfdim, in, out, FFTW_FORWARD, FFTW_MEASURE);


    /*Number of grid points in th enon-local box in coarse and double grids */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    weptr1 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr2 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr3 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    gwptr = (double complex *)fftw_malloc(sizeof(double complex) * size);

    if ((weptr1 == NULL) || (weptr2 == NULL) || (weptr3 == NULL) || (gwptr == NULL))
        error_handler ("can't allocate memory\n");

    hxx = get_hxgrid() / (double) ct.nxfgrid;
    hyy = get_hygrid() / (double) ct.nyfgrid;
    hzz = get_hzgrid() / (double) ct.nzfgrid;

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
                //t1 = linint (&sp->betalig[ip][0], r, invdr);
                t1 = AtomicInterpolate (&sp->betalig[ip][0], r);
                to_cartesian (ax, bx);
                r += 1.0e-10;

                weptr1[idx] = cc * bx[0] * t1 / r + 0.0I;
                weptr2[idx] = cc * bx[2] * t1 / r + 0.0I;
                weptr3[idx] = cc * bx[1] * t1 / r + 0.0I;


                if((ix*2 + sp->nlfdim) == 0 || (iy*2 + sp->nlfdim) == 0 || (iz*2 + sp->nlfdim) == 0 )
                {
                    weptr1[idx] = 0.0;
                    weptr2[idx] = 0.0;
                    weptr3[idx] = 0.0;
                }


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    int broot[3], jdx;
    int npes = get_PE_X() * get_PE_Y() * get_PE_Z();
    int istop = 3;
    if(npes < istop) istop = npes;
    idx = 0;
    while(idx < 3) {
        for(jdx = 0;jdx < istop;jdx++) {
            broot[idx] = jdx;
            idx++;
        }
    }

    if(pct.gridpe == broot[0]) {
        fftw_execute_dft (p1, weptr1, gwptr);
        pack_gftoc (sp, gwptr, r1);
    }

    if(pct.gridpe == broot[1]) {
        fftw_execute_dft (p1, weptr2, gwptr);
        pack_gftoc (sp, gwptr, r2);
    }

    if(pct.gridpe == broot[2]) {
        fftw_execute_dft (p1, weptr3, gwptr);
        pack_gftoc (sp, gwptr, r3);
    }

    MPI_Bcast(r1, 2*coarse_size, MPI_DOUBLE, broot[0], pct.grid_comm);
    MPI_Bcast(r2, 2*coarse_size, MPI_DOUBLE, broot[1], pct.grid_comm);
    MPI_Bcast(r3, 2*coarse_size, MPI_DOUBLE, broot[2], pct.grid_comm);


    fftw_free (gwptr);
    fftw_free (weptr3);
    fftw_free (weptr2);
    fftw_free (weptr1);

}
