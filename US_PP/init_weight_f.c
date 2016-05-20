/************************** SVN Revision Information **************************
 **    $Id: init_weight_d.c 3588 2016-05-20 12:29:40Z ebriggs $    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "main.h"
#include "common_prototypes.h"
#include "AtomicInterpolate.h"

void init_weight_f (SPECIES * sp, fftw_complex * rtptr, int ip, fftw_plan p1, bool use_shared)
{

    if(use_shared && (pct.local_rank > 6)) return;

    int idx, ix, iy, iz, size, coarse_size, iend, ibegin;
    double r, ax[3], bx[3], xc, yc, zc, t1, t2, rsq1;
    double c0, c1, c2, c3, hxx, hyy, hzz;
    double complex *weptr1, *weptr2, *weptr3, *weptr4, *weptr5, *weptr6, *weptr7, *gwptr;
    double complex *r1, *r2, *r3, *r4, *r5, *r6, *r7;
    int ixx, iyy, izz;


    /*Number of grid points in th enon-local box in coarse and double grids */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    weptr1 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr2 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr3 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr4 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr5 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr6 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr7 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    gwptr = (double complex *)fftw_malloc(sizeof(double complex) * size);

    if ((weptr1 == NULL) || (weptr2 == NULL) || (weptr3 == NULL) || (weptr4 == NULL) ||
        (weptr5 == NULL) || (weptr6 == NULL) || (weptr7 == NULL) || (gwptr == NULL))
        error_handler ("can't allocate memory\n");

    hxx = get_hxgrid() / (double) ct.nxfgrid;
    hyy = get_hygrid() / (double) ct.nyfgrid;
    hzz = get_hzgrid() / (double) ct.nzfgrid;

    r1 = rtptr;
    r2 = r1 + coarse_size;
    r3 = r2 + coarse_size;
    r4 = r3 + coarse_size;
    r5 = r4 + coarse_size;
    r6 = r5 + coarse_size;
    r7 = r6 + coarse_size;

    c0 = sqrt (35.0 / (2.0 * PI)) / 4.0;
    c1 = sqrt(105.0 / PI);
    c2 = sqrt(21.0 / (2.0 * PI)) / 4.0;
    c3 = sqrt(7.0 / PI) / 4.0;

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

                rsq1 = r * r * r + 1.0e-20;
                t1 = AtomicInterpolateInline (&sp->betalig[ip][0], r);
                weptr1[idx] = c0 * t1 * (bx[1] * (3.0*bx[0]*bx[0] - bx[1]*bx[1])) / rsq1 + 0.0I;
                weptr2[idx] = c1 * t1 * (bx[0] * bx[1] * bx[2]) / (2.0*rsq1) + 0.0I;
                weptr3[idx] = c2 * t1 * (bx[1] * (4.0*bx[2]*bx[2] - bx[0]*bx[0] - bx[1]*bx[1])) / rsq1 + 0.0I;
                weptr4[idx] = c3 * t1 * (bx[2] * (2.0*bx[2]*bx[2] - 3.0*bx[0]*bx[0] - 3.0*bx[1]*bx[1])) / rsq1 + 0.0I;
                weptr5[idx] = c2 * t1 * (bx[0] * (4.0*bx[2]*bx[2] - bx[0]*bx[0] - bx[1]*bx[1])) / rsq1 + 0.0I;
                weptr6[idx] = c1 * t1 * (bx[2] * (bx[0]*bx[0] - bx[1]*bx[1])) / (4.0*rsq1) + 0.0I;
                weptr7[idx] = c0 * t1 * (bx[0] * (bx[0]*bx[0] - 3.0*bx[1]*bx[1])) / rsq1 + 0.0I;


                if((ix*2 + sp->nlfdim) == 0 || (iy*2 + sp->nlfdim) == 0 || (iz*2 + sp->nlfdim) == 0 )
                {
                    weptr1[idx] = 0.0;
                    weptr2[idx] = 0.0;
                    weptr3[idx] = 0.0;
                    weptr4[idx] = 0.0;
                    weptr5[idx] = 0.0;
                    weptr6[idx] = 0.0;
                    weptr7[idx] = 0.0;
                }

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    if(use_shared) {

        // Always at least 7 procs per host if use_shared is true
        if(pct.local_rank == 0) {
            fftw_execute_dft (p1, weptr1, gwptr);
            pack_gftoc (sp, gwptr, r1);
        }
        if(pct.local_rank == 1) {
            fftw_execute_dft (p1, weptr2, gwptr);
            pack_gftoc (sp, gwptr, r2);
        }
        if(pct.local_rank == 2) {
            fftw_execute_dft (p1, weptr3, gwptr);
            pack_gftoc (sp, gwptr, r3);
        }
        if(pct.local_rank == 3) {
            fftw_execute_dft (p1, weptr4, gwptr);
            pack_gftoc (sp, gwptr, r4);
        }
        if(pct.local_rank == 4) {
            fftw_execute_dft (p1, weptr5, gwptr);
            pack_gftoc (sp, gwptr, r5);
        }
        if(pct.local_rank == 5) {
            fftw_execute_dft (p1, weptr6, gwptr);
            pack_gftoc (sp, gwptr, r6);
        }
        if(pct.local_rank == 6) {
            fftw_execute_dft (p1, weptr7, gwptr);
            pack_gftoc (sp, gwptr, r7);
        }

    }
    else {

        int broot[7], jdx;
        int npes = get_PE_X() * get_PE_Y() * get_PE_Z();
        int istop = 7;
        if(npes < istop) istop = npes;
        for(idx = 0; idx < 7; idx++)
            broot[idx] = idx%istop;

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

        if(pct.gridpe == broot[3]) {
            fftw_execute_dft (p1, weptr4, gwptr);
            pack_gftoc (sp, gwptr, r4);
        }

        if(pct.gridpe == broot[4]) {
            fftw_execute_dft (p1, weptr5, gwptr);
            pack_gftoc (sp, gwptr, r5);
        }

        if(pct.gridpe == broot[5]) {
            fftw_execute_dft (p1, weptr6, gwptr);
            pack_gftoc (sp, gwptr, r6);
        }

        if(pct.gridpe == broot[6]) {
            fftw_execute_dft (p1, weptr7, gwptr);
            pack_gftoc (sp, gwptr, r7);
        }

        MPI_Bcast(r1, 2*coarse_size, MPI_DOUBLE, broot[0], pct.grid_comm);
        MPI_Bcast(r2, 2*coarse_size, MPI_DOUBLE, broot[1], pct.grid_comm);
        MPI_Bcast(r3, 2*coarse_size, MPI_DOUBLE, broot[2], pct.grid_comm);
        MPI_Bcast(r4, 2*coarse_size, MPI_DOUBLE, broot[3], pct.grid_comm);
        MPI_Bcast(r5, 2*coarse_size, MPI_DOUBLE, broot[4], pct.grid_comm);
        MPI_Bcast(r6, 2*coarse_size, MPI_DOUBLE, broot[5], pct.grid_comm);
        MPI_Bcast(r7, 2*coarse_size, MPI_DOUBLE, broot[6], pct.grid_comm);

    }


    fftw_free (gwptr);
    fftw_free (weptr7);
    fftw_free (weptr6);
    fftw_free (weptr5);
    fftw_free (weptr4);
    fftw_free (weptr3);
    fftw_free (weptr2);
    fftw_free (weptr1);

}
