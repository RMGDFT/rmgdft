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

void init_weight_s (SPECIES * sp, fftw_complex * rtptr, int ip, fftw_plan p1, bool use_shared)
{

    if(use_shared && (pct.local_rank != 0)) return;
    int idx, ix, iy, iz, size;
    double r, ax[3], xc, yc, zc;
    double t1;
    double complex *weptr, *gwptr;

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

    weptr = (double complex *)fftw_malloc(sizeof(double complex) * size);
    gwptr = (double complex *)fftw_malloc(sizeof(double complex) * size);

    if ((weptr == NULL) || (gwptr == NULL))
        error_handler ("can't allocate memory\n");


    /*We assume that ion is in the center of non-local box */
//    ibegin = -(sp->nldim / 2) * ct.nxfgrid;


    int ixbegin = -nlfxdim/2;
    int ixend = ixbegin + nlfxdim;
    int iybegin = -nlfydim/2;
    int iyend = iybegin + nlfydim;
    int izbegin = -nlfzdim/2;
    int izend = izbegin + nlfzdim;

    int ixx, iyy, izz;
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

                ax[0] = xc;
                ax[1] = yc;
                ax[2] = zc;

                r = metric (ax);

                t1 = AtomicInterpolateInline(&sp->betalig[ip][0], r);
                idx = ixx * nlfydim * nlfzdim + iyy * nlfzdim + izz;
                weptr[idx] = sqrt (1.0 / (4.0 * PI)) * t1 + 0.0I;

                if((ix*2 + sp->nlfdim) == 0 || (iy*2 + sp->nlfdim) == 0 || (iz*2 + sp->nlfdim) == 0 ) 
                    weptr[idx] = 0.0;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    fftw_execute_dft (p1, weptr, gwptr);

    pack_gftoc (sp, gwptr, rtptr);

    fftw_free (gwptr);
    fftw_free (weptr);
}
