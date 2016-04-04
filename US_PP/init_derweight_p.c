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
#include "AtomicInterpolate.h"

void init_derweight_p (SPECIES * sp,
                       fftw_complex * rtptr_x,
                       fftw_complex * rtptr_y, fftw_complex * rtptr_z, int ip, fftw_plan p1)
{

    int idx, ix, iy, iz, size, coarse_size, ibegin, iend;
    double r, rsq, r3, rsqd, ax[3], bx[3], x, y, z, xc, yc, zc, cc, t1, t2;
    double hxx, hyy, hzz;
    double complex *weptr1x, *weptr1y, *weptr1z, *gwptr;
    double complex *weptr2x, *weptr2y, *weptr2z;
    double complex *weptr3x, *weptr3y, *weptr3z;
    double complex *r1x, *r1y, *r1z, *r2x, *r2y, *r2z, *r3x, *r3y, *r3z;
    int ixx, iyy, izz;

    /*Number of grid points in the non-local box in coarse and double grids */
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

    hxx = get_hxgrid() / (double) ct.nxfgrid;
    hyy = get_hygrid() / (double) ct.nyfgrid;
    hzz = get_hzgrid() / (double) ct.nzfgrid;

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

                t1 = AtomicInterpolateInline (&sp->drbetalig[ip][0], r);
                t2 = AtomicInterpolateInline (&sp->betalig[ip][0], r);

                to_cartesian (ax, bx);
                x = bx[0];
                y = bx[1];
                z = bx[2];

                rsq = x * x + y * y + z * z;
                r3 = r * rsq + 1.0e-20;
                rsqd = rsq + 1.0e-20;


                weptr1x[idx] = cc * ((t2 * (rsq - x * x) / r3) + (x * x * t1 / rsqd)) + 0.0I;
                weptr1y[idx] = cc * ((-t2 * x * y / r3) + (x * y * t1 / rsqd)) + 0.0I;
                weptr1z[idx] = cc * ((-t2 * x * z / r3) + (x * z * t1 / rsqd)) + 0.0I;

                weptr2x[idx] = cc * ((-t2 * x * z / r3) + (x * z * t1 / rsqd)) + 0.0I;
                weptr2y[idx] = cc * ((-t2 * y * z / r3) + (y * z * t1 / rsqd)) + 0.0I;
                weptr2z[idx] = cc * ((t2 * (rsq - z * z) / r3) + (z * z * t1 / rsqd)) + 0.0I;

                weptr3x[idx] = cc * ((-t2 * x * y / r3) + (x * y * t1 / rsqd)) + 0.0I;
                weptr3y[idx] = cc * ((t2 * (rsq - y * y) / r3) + (y * y * t1 / rsqd)) + 0.0I;
                weptr3z[idx] = cc * ((-t2 * y * z / r3) + (y * z * t1 / rsqd)) + 0.0I;



                if( ((ix*2 + sp->nlfdim) == 0) || ((iy*2 + sp->nlfdim) == 0) || ((iz*2 + sp->nlfdim) == 0) )
                {
                    weptr1x[idx] = 0.0;
                    weptr1y[idx] = 0.0;
                    weptr1z[idx] = 0.0;
                    weptr2x[idx] = 0.0;
                    weptr2y[idx] = 0.0;
                    weptr2z[idx] = 0.0;
                    weptr3x[idx] = 0.0;
                    weptr3y[idx] = 0.0;
                    weptr3z[idx] = 0.0;
                }

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    int broot[9], jdx;
    int npes = get_PE_X() * get_PE_Y() * get_PE_Z();
    int istop = 9;
    if(npes < istop) istop = npes;
    for(idx = 0; idx < 9; idx++)
        broot[idx] = idx %istop;

    if(pct.gridpe == broot[0]) {
        fftw_execute_dft (p1, weptr1x, gwptr);
        pack_gftoc (sp, gwptr, r1x);
    }

    if(pct.gridpe == broot[1]) {
        fftw_execute_dft (p1, weptr1y, gwptr);
        pack_gftoc (sp, gwptr, r1y);
    }

    if(pct.gridpe == broot[2]) {
        fftw_execute_dft (p1, weptr1z, gwptr);
        pack_gftoc (sp, gwptr, r1z);
    }


    if(pct.gridpe == broot[3]) {
        fftw_execute_dft (p1, weptr2x, gwptr);
        pack_gftoc (sp, gwptr, r2x);
    }

    if(pct.gridpe == broot[4]) {
        fftw_execute_dft (p1, weptr2y, gwptr);
        pack_gftoc (sp, gwptr, r2y);
    }

    if(pct.gridpe == broot[5]) {
        fftw_execute_dft (p1, weptr2z, gwptr);
        pack_gftoc (sp, gwptr, r2z);
    }


    if(pct.gridpe == broot[6]) {
        fftw_execute_dft (p1, weptr3x, gwptr);
        pack_gftoc (sp, gwptr, r3x);
    }

    if(pct.gridpe == broot[7]) {
        fftw_execute_dft (p1, weptr3y, gwptr);
        pack_gftoc (sp, gwptr, r3y);
    }

    if(pct.gridpe == broot[8]) {
        fftw_execute_dft (p1, weptr3z, gwptr);
        pack_gftoc (sp, gwptr, r3z);
    }


    MPI_Bcast(r1x, 2*coarse_size, MPI_DOUBLE, broot[0], pct.grid_comm);
    MPI_Bcast(r1y, 2*coarse_size, MPI_DOUBLE, broot[1], pct.grid_comm);
    MPI_Bcast(r1z, 2*coarse_size, MPI_DOUBLE, broot[2], pct.grid_comm);
    MPI_Bcast(r2x, 2*coarse_size, MPI_DOUBLE, broot[3], pct.grid_comm);
    MPI_Bcast(r2y, 2*coarse_size, MPI_DOUBLE, broot[4], pct.grid_comm);
    MPI_Bcast(r2z, 2*coarse_size, MPI_DOUBLE, broot[5], pct.grid_comm);
    MPI_Bcast(r3x, 2*coarse_size, MPI_DOUBLE, broot[6], pct.grid_comm);
    MPI_Bcast(r3y, 2*coarse_size, MPI_DOUBLE, broot[7], pct.grid_comm);
    MPI_Bcast(r3z, 2*coarse_size, MPI_DOUBLE, broot[8], pct.grid_comm);



    my_free (weptr1x);

}
