/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
                            make_mask_grid.c

    make mask_grid for each state according their radius and mask level  

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void make_mask_grid_state(int level, STATE * states)
{
    int ix, iy, iz;
    int ii;
    REAL r, xc[3], offset[3], hgrid[3], delta;
    char *ptr_mask, *maskptr;
    int dim[3], incx, incy;
    char myMask;
    int state;
    REAL rcut;


    hgrid[0] = ct.hxgrid * ct.xside * (REAL) (1 << level);
    hgrid[1] = ct.hygrid * ct.yside * (REAL) (1 << level);
    hgrid[2] = ct.hzgrid * ct.zside * (REAL) (1 << level);

    /* Coeff. 2 for the 4th order FD scheme in postprocessing */
    /*delta = 2.*sqrt( hgrid[0]*hgrid[0]+hgrid[1]*hgrid[1]+hgrid[2]*hgrid[2]); */
    delta = 2. * hgrid[0];
    delta = max(delta, hgrid[1]);
    delta = max(delta, hgrid[2]);

    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        rcut = states[state].radius;
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
/*		if(myMask == 1) printf("\n state, x, y, z: %d, %d, %d, %d,  mask : %d", state, ix, iy, iz, myMask); 
 *		*/
                }
                xc[1] += hgrid[1];
            }
            xc[0] += hgrid[0];
        }
    }

    if (pct.gridpe == 0)
        printf(" make_mask_grid_state is  done\n");
}
