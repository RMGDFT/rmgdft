/************************** SVN Revision Information **************************
 **    $Id: partial_vlocpsi.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/*
partial_vlocpsi.c

multiplication of (orbit,  local projector )
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void partial_vlocpsi (STATE st1, int ion2, REAL * psi, REAL * prjptr, REAL *vlpsi)
{

    int xlow1, xhigh1, xlow2, xhigh2, xshift;
    int ylow1, yhigh1, ylow2, yhigh2, yshift;
    int zlow1, zhigh1, zlow2, zhigh2, zshift;
    int iyy, izz, iyy1, izz1;
    int incx, incy, incx1, incy1;
    int ix, iy, iz, ix1, ix2, iy1, iy2, iz1, iz2;
    int idx1, idx2;
    int index;

    REAL time1, time2;
    time1 = my_crtc ();

/*    for (idx1=0; idx1 < st1.size; idx1++) vlpsi[idx1] = 0.0;*/

    index = (st1.index - ct.state_begin) * ct.num_ions + ion2;

    xlow1 = ion_orbit_overlap_region_loc[index].xlow1;
    xhigh1 = ion_orbit_overlap_region_loc[index].xhigh1;
    xlow2 = ion_orbit_overlap_region_loc[index].xlow2;
    xhigh2 = ion_orbit_overlap_region_loc[index].xhigh2;
    xshift = ion_orbit_overlap_region_loc[index].xshift;

    ylow1 = ion_orbit_overlap_region_loc[index].ylow1;
    yhigh1 = ion_orbit_overlap_region_loc[index].yhigh1;
    ylow2 = ion_orbit_overlap_region_loc[index].ylow2;
    yhigh2 = ion_orbit_overlap_region_loc[index].yhigh2;
    yshift = ion_orbit_overlap_region_loc[index].yshift;

    zlow1 = ion_orbit_overlap_region_loc[index].zlow1;
    zhigh1 = ion_orbit_overlap_region_loc[index].zhigh1;
    zlow2 = ion_orbit_overlap_region_loc[index].zlow2;
    zhigh2 = ion_orbit_overlap_region_loc[index].zhigh2;
    zshift = ion_orbit_overlap_region_loc[index].zshift;

    iyy = st1.iymax - st1.iymin + 1;
    izz = st1.izmax - st1.izmin + 1;
    incx = iyy * izz;
    incy = izz;

    iyy1 = ct.ions[ion2].iyend - ct.ions[ion2].iystart + 1;
    izz1 = ct.ions[ion2].izend - ct.ions[ion2].izstart + 1;
    incx1 = iyy1 * izz1;
    incy1 = izz1;


    for (ix = xlow1; ix <= xhigh1; ix++)
    {
        ix1 = (ix - st1.ixmin) * incx;
        ix2 = (ix - ct.ions[ion2].ixstart) * incx1;

        for (iy = ylow1; iy <= yhigh1; iy++)
        {
            iy1 = (iy - st1.iymin) * incy;
            iy2 = (iy - ct.ions[ion2].iystart) * incy1;

            for (iz = zlow1; iz <= zhigh1; iz++)
            {
                iz1 = iz - st1.izmin;
                iz2 = iz - ct.ions[ion2].izstart;
                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                vlpsi[idx1] = psi[idx1] * prjptr[idx2];
            }

            for (iz = zlow2; iz <= zhigh2; iz++)
            {
                iz1 = iz - st1.izmin;
                iz2 = iz - ct.ions[ion2].izstart - zshift;
                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                vlpsi[idx1] = psi[idx1] * prjptr[idx2];
            }
        }

        for (iy = ylow2; iy <= yhigh2; iy++)
        {
            iy1 = (iy - st1.iymin) * incy;
            iy2 = (iy - ct.ions[ion2].iystart - yshift) * incy1;

            for (iz = zlow1; iz <= zhigh1; iz++)
            {
                iz1 = iz - st1.izmin;
                iz2 = iz - ct.ions[ion2].izstart;
                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                vlpsi[idx1] = psi[idx1] * prjptr[idx2];
            }
            for (iz = zlow2; iz <= zhigh2; iz++)
            {
                iz1 = iz - st1.izmin;
                iz2 = iz - ct.ions[ion2].izstart - zshift;
                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                vlpsi[idx1] = psi[idx1] * prjptr[idx2];
            }
        }
    }                           /* end for ix = xlow1 */

    for (ix = xlow2; ix <= xhigh2; ix++)
    {
        ix1 = (ix - st1.ixmin) * incx;
        ix2 = (ix - ct.ions[ion2].ixstart - xshift) * incx1;

        for (iy = ylow1; iy <= yhigh1; iy++)
        {
            iy1 = (iy - st1.iymin) * incy;
            iy2 = (iy - ct.ions[ion2].iystart) * incy1;
            for (iz = zlow1; iz <= zhigh1; iz++)
            {
                iz1 = iz - st1.izmin;
                iz2 = iz - ct.ions[ion2].izstart;
                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                vlpsi[idx1] = psi[idx1] * prjptr[idx2];
            }
            for (iz = zlow2; iz <= zhigh2; iz++)
            {
                iz1 = iz - st1.izmin;
                iz2 = iz - ct.ions[ion2].izstart - zshift;
                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                vlpsi[idx1] = psi[idx1] * prjptr[idx2];
            }
        }                       /* end for iy = ylow1 */

        for (iy = ylow2; iy <= yhigh2; iy++)
        {
            iy1 = (iy - st1.iymin) * incy;
            iy2 = (iy - ct.ions[ion2].iystart - yshift) * incy1;
            for (iz = zlow1; iz <= zhigh1; iz++)
            {
                iz1 = iz - st1.izmin;
                iz2 = iz - ct.ions[ion2].izstart;
                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                vlpsi[idx1] = psi[idx1] * prjptr[idx2];
            }
            for (iz = zlow2; iz <= zhigh2; iz++)
            {
                iz1 = iz - st1.izmin;
                iz2 = iz - ct.ions[ion2].izstart - zshift;
                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                vlpsi[idx1] = psi[idx1] * prjptr[idx2];
            }
        }
    }                           /* end for ix = xlow2 */

    time2 = my_crtc ();
    md_timings (PAR_VLOCPSI, (time2 - time1));

}
