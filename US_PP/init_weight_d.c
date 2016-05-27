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
#include "AtomicInterpolate.h"

void init_weight_d (SPECIES * sp, fftw_complex * rtptr, int ip, fftw_plan p1, bool use_shared)
{

    if(use_shared && (pct.local_rank > 4)) return;

    int idx, ix, iy, iz, size, coarse_size;
    double r, ax[3], bx[3], xc, yc, zc, t1, t2, rsq1;
    double cc;
    double complex *weptr1, *weptr2, *weptr3, *weptr4, *weptr5, *gwptr;
    double complex *r1, *r2, *r3, *r4, *r5;
    int ixx, iyy, izz;

    double hxx = get_hxgrid() / (double) ct.nxfgrid;
    double hyy = get_hygrid() / (double) ct.nyfgrid;
    double hzz = get_hzgrid() / (double) ct.nzfgrid;
    double xoff = 0.0;
    double yoff = 0.0;
    double zoff = 0.0;

    int nlfxdim = sp->nlfdim;
    int nlfydim = sp->nlfdim;
    int nlfzdim = sp->nlfdim;
    if(!ct.localize_projectors) {
        nlfxdim = get_NX_GRID();
        nlfydim = get_NY_GRID();
        nlfzdim = get_NZ_GRID();
        xoff = 0.5 * hxx;
        yoff = 0.5 * hyy;
        zoff = 0.5 * hzz;
    }

    /*Number of grid points in the non-local box in coarse and double grids */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;
    size = nlfxdim * nlfydim * nlfzdim;

    weptr1 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr2 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr3 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr4 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr5 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    gwptr = (double complex *)fftw_malloc(sizeof(double complex) * size);

    if ((weptr1 == NULL) || (weptr2 == NULL) || (weptr3 == NULL) || (weptr4 == NULL) || (weptr5 == NULL) || (gwptr == NULL))
        error_handler ("can't allocate memory\n");

    r1 = rtptr;
    r2 = r1 + coarse_size;
    r3 = r2 + coarse_size;
    r4 = r3 + coarse_size;
    r5 = r4 + coarse_size;

    cc = sqrt (5.0 / (4.0 * PI));
    t2 = sqrt (3.0);


    int ixbegin = -nlfxdim/2;
    int ixend = ixbegin + nlfxdim;
    int iybegin = -nlfydim/2;
    int iyend = iybegin + nlfydim;
    int izbegin = -nlfzdim/2;
    int izend = izbegin + nlfzdim;

    for (ix = ixbegin; ix < ixend; ix++)
    {
        ixx = ix;
        if (ixx < 0) ixx = ix + nlfxdim;
        xc = (double) ix *hxx + xoff;

        for (iy = iybegin; iy < iyend; iy++)
        {
            iyy = iy;
            if (iyy < 0) iyy = iy + nlfydim;
            yc = (double) iy *hyy + yoff;

            for (iz = izbegin; iz < izend; iz++)
            {

                izz = iz;
                if (izz < 0) izz = iz + nlfzdim;
                zc = (double) iz *hzz + zoff;

                idx = ixx * nlfxdim * nlfydim + iyy * nlfzdim + izz;

                ax[0] = xc;
                ax[1] = yc;
                ax[2] = zc;

                to_cartesian (ax, bx);
                r = metric (ax);

                rsq1 = r * r + 1.0e-20;
                t1 = AtomicInterpolateInline (&sp->betalig[ip][0], r);
                t1 = t1 * t2;
                weptr1[idx] = cc * t1 * bx[0] * bx[1] / rsq1 + 0.0I;
                weptr2[idx] = cc * t1 * bx[0] * bx[2] / rsq1 + 0.0I;
                weptr3[idx] = cc * t1 * (t2 * bx[2] * bx[2] - rsq1 / t2) / (2.0 * rsq1) + 0.0I;
                weptr4[idx] = cc * t1 * bx[1] * bx[2] / rsq1 + 0.0I;
                weptr5[idx] = cc * t1 * (bx[0] * bx[0] - bx[1] * bx[1]) / (2.0 * rsq1) + 0.0I;


                if((ix*2 + sp->nlfdim) == 0 || (iy*2 + sp->nlfdim) == 0 || (iz*2 + sp->nlfdim) == 0 )
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

    if(use_shared) {

        // Always at least 5 procs per host if use_shared is true
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

    }
    else {

        int broot[5];
        int npes = get_PE_X() * get_PE_Y() * get_PE_Z();
        int istop = 5;
        if(npes < istop) istop = npes;
        for(idx = 0; idx < 5; idx++)
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

        MPI_Bcast(r1, 2*coarse_size, MPI_DOUBLE, broot[0], pct.grid_comm);
        MPI_Bcast(r2, 2*coarse_size, MPI_DOUBLE, broot[1], pct.grid_comm);
        MPI_Bcast(r3, 2*coarse_size, MPI_DOUBLE, broot[2], pct.grid_comm);
        MPI_Bcast(r4, 2*coarse_size, MPI_DOUBLE, broot[3], pct.grid_comm);
        MPI_Bcast(r5, 2*coarse_size, MPI_DOUBLE, broot[4], pct.grid_comm);

    }


    fftw_free (gwptr);
    fftw_free (weptr5);
    fftw_free (weptr4);
    fftw_free (weptr3);
    fftw_free (weptr2);
    fftw_free (weptr1);

}
