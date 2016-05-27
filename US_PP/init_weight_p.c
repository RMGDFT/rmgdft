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


void init_weight_p (SPECIES * sp, fftw_complex * rtptr, int ip, fftw_plan p1, bool use_shared)
{

    if(use_shared && (pct.local_rank > 2)) return;

    int idx, ix, iy, iz, size, coarse_size;
    double r, ax[3], bx[3], xc, yc, zc, cc, t1;
    double complex *weptr1, *weptr2, *weptr3, *gwptr;
    double complex *r1, *r2, *r3;
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
        nlfxdim = ct.nxfgrid * get_NX_GRID();
        nlfydim = ct.nxfgrid * get_NY_GRID();
        nlfzdim = ct.nxfgrid * get_NZ_GRID();
        xoff = 0.5 * hxx;
        yoff = 0.5 * hyy;
        zoff = 0.5 * hzz;
    }

    /* nlxdim * nlydim * nlzdim is size of the non-local box in the double grid */
    size = nlfxdim * nlfydim * nlfzdim;

    /*Number of grid points in th enon-local box in coarse and double grids */
    coarse_size = sp->nldim * sp->nldim * sp->nldim;

    weptr1 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr2 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    weptr3 = (double complex *)fftw_malloc(sizeof(double complex) * size);
    gwptr = (double complex *)fftw_malloc(sizeof(double complex) * size);

    if ((weptr1 == NULL) || (weptr2 == NULL) || (weptr3 == NULL) || (gwptr == NULL))
        error_handler ("can't allocate memory\n");

    r1 = rtptr;
    r2 = r1 + coarse_size;
    r3 = r2 + coarse_size;

    cc = sqrt (3.0 / (4.0 * PI));

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

                idx = ixx *sp->nlfdim * sp->nlfdim + iyy * sp->nlfdim + izz;

                ax[0] = xc;
                ax[1] = yc;
                ax[2] = zc;

                r = metric (ax);
                t1 = AtomicInterpolateInline (&sp->betalig[ip][0], r);
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

    if(use_shared) {

        // Always at least 3 procs per host if use_shared is true
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

    }
    else {

        int broot[3];
        int npes = get_PE_X() * get_PE_Y() * get_PE_Z();
        int istop = 3;
        if(npes < istop) istop = npes;
        for(idx = 0; idx < 3; idx++)
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

        MPI_Bcast(r1, 2*coarse_size, MPI_DOUBLE, broot[0], pct.grid_comm);
        MPI_Bcast(r2, 2*coarse_size, MPI_DOUBLE, broot[1], pct.grid_comm);
        MPI_Bcast(r3, 2*coarse_size, MPI_DOUBLE, broot[2], pct.grid_comm);

    }
    

    fftw_free (gwptr);
    fftw_free (weptr3);
    fftw_free (weptr2);
    fftw_free (weptr1);

}
