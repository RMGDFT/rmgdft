/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
                            make_mask_grid.c
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void make_mask_grid(REAL rcut, int level, STATE * states)
{
    int ix, iy, iz;
    int ii;
    REAL r, xc[3], offset[3], hgrid[3], delta;
    char *ptr_mask, *maskptr;
    int dim[3], incx, incy;
    char myMask;
    int state;


    hgrid[0] = ct.hxgrid * ct.xside * (REAL) (1 << level);
    hgrid[1] = ct.hygrid * ct.yside * (REAL) (1 << level);
    hgrid[2] = ct.hzgrid * ct.zside * (REAL) (1 << level);

    /* Coeff. 2 for the 4th order FD scheme in postprocessing */
    /*delta = 2.*sqrt( hgrid[0]*hgrid[0]+hgrid[1]*hgrid[1]+hgrid[2]*hgrid[2]); */
    delta = 2. * hgrid[0];
    delta = max(delta, hgrid[1]);
    delta = max(delta, hgrid[2]);

    if (pct.gridpe == 0)
    {
        printf("\n  make_mask_grid: level %d, rcut = %f rcut+delta = %f\n",
               level, rcut, rcut + delta);
        if (ct.boundaryflag == ZPERIODIC)
            printf(" make_mask_grid: ZPERIODIC\n");
    }

    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        dim[0] = states[state].orbit_nx / (1 << level);
        dim[1] = states[state].orbit_ny / (1 << level);
        dim[2] = states[state].orbit_nz / (1 << level);
        incy = dim[2];
        incx = dim[2] * dim[1];

        ptr_mask = states[state].lmask[0];
        for (ii = 0; ii < level; ii++)
            ptr_mask += states[state].size / ((1 << ii) * (1 << ii) * (1 << ii));
        states[state].lmask[level] = ptr_mask;
        maskptr = ptr_mask;

        offset[0] = states[state].ixmin * ct.hxgrid * ct.xside;
        offset[1] = states[state].iymin * ct.hygrid * ct.yside;
        offset[2] = states[state].izmin * ct.hzgrid * ct.zside;

        xc[0] = offset[0];
        for (ix = 0; ix < dim[0]; ix++)
        {
            xc[1] = offset[1];
            for (iy = 0; iy < dim[1]; iy++)
            {
                xc[2] = offset[2];
                for (iz = 0; iz < dim[2]; iz++)
                {
                    r = minimage1(xc, states[state].crds);
                    if (r > rcut)
                        myMask = 0;
                    else
                        myMask = 1;
                    maskptr[ix * incx + iy * incy + iz] = myMask;
                    xc[2] += hgrid[2];
                }
                xc[1] += hgrid[1];
            }
            xc[0] += hgrid[0];
        }
    }

    if (pct.gridpe == 0)
        printf(" make_mask_grid: mask done\n");
}



REAL minimage1(REAL aa[3], REAL bb[3])
{

    int ix, iy, iz, idx, idxmin, nn = 0;
    REAL r[27], ax, ay, az, rmin, x, y, z;


    /* Get lattice vectors */
    ax = NX_GRID * ct.hxgrid * ct.xside;
    ay = NY_GRID * ct.hygrid * ct.yside;
    az = NZ_GRID * ct.hzgrid * ct.zside;



    /* Loop over all possible combinations and generate r */
    if (ct.boundaryflag == PERIODIC)
    {
        nn = 27;

        idx = 0;
        for (ix = -1; ix <= 1; ix++)
        {

            for (iy = -1; iy <= 1; iy++)
            {

                for (iz = -1; iz <= 1; iz++)
                {

                    x = aa[0] - bb[0] + (REAL) ix *ax;
                    y = aa[1] - bb[1] + (REAL) iy *ay;
                    z = aa[2] - bb[2] + (REAL) iz *az;

                    r[idx] = (x * x + y * y + z * z);

                    idx++;

                }               /* end for */

            }                   /* end for */

        }                       /* end for */
    }
    else if (ct.boundaryflag == ZPERIODIC)
    {
        nn = 3;

        idx = 0;
        for (iz = -1; iz <= 1; iz++)
        {

            x = aa[0] - bb[0];
            y = aa[1] - bb[1];
            z = aa[2] - bb[2] + (REAL) iz *az;

            r[idx] = (x * x + y * y + z * z);

            idx++;

        }                       /* end for */
    }


    /* Next we find rmin */
    rmin = 10000000.0;
    for (idx = 0; idx < nn; idx++)
    {

        if (r[idx] < rmin)
        {

            rmin = r[idx];
            idxmin = idx;

        }                       /* end if */

    }                           /* end for */

    rmin = sqrt(rmin);
    return rmin;

}                               /* end minimage1 */
