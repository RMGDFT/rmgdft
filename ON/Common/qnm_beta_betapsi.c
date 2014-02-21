/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*


   kb_kb_orbit.c

   calculate |kb><kb|psi>

 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main_on.h"

/*shuchun wang */
void qnm_beta_betapsi(STATE *state1, int ion2, rmg_double_t * pptr, rmg_double_t * work)
{

    int xlow1, xhigh1, xlow2, xhigh2, xshift;
    int ylow1, yhigh1, ylow2, yhigh2, yshift;
    int zlow1, zhigh1, zlow2, zhigh2, zshift;
    int iyy, izz, iyy1, izz1;
    int incx, incy, incx1, incy1;
    int ix, iy, iz;
    int ix1, ix2, iy1, iy2, iz1, iz2, idx1, idx2;
    int index;


    index = (state1->index - ct.state_begin) * ct.num_ions + ion2;

    xlow1 = ion_orbit_overlap_region_nl[index].xlow1;
    xhigh1 = ion_orbit_overlap_region_nl[index].xhigh1;
    xlow2 = ion_orbit_overlap_region_nl[index].xlow2;
    xhigh2 = ion_orbit_overlap_region_nl[index].xhigh2;
    xshift = ion_orbit_overlap_region_nl[index].xshift;

    ylow1 = ion_orbit_overlap_region_nl[index].ylow1;
    yhigh1 = ion_orbit_overlap_region_nl[index].yhigh1;
    ylow2 = ion_orbit_overlap_region_nl[index].ylow2;
    yhigh2 = ion_orbit_overlap_region_nl[index].yhigh2;
    yshift = ion_orbit_overlap_region_nl[index].yshift;

    zlow1 = ion_orbit_overlap_region_nl[index].zlow1;
    zhigh1 = ion_orbit_overlap_region_nl[index].zhigh1;
    zlow2 = ion_orbit_overlap_region_nl[index].zlow2;
    zhigh2 = ion_orbit_overlap_region_nl[index].zhigh2;
    zshift = ion_orbit_overlap_region_nl[index].zshift;

    iyy = state1->iymax - state1->iymin + 1;
    izz = state1->izmax - state1->izmin + 1;
    incx = iyy * izz;
    incy = izz;

    iyy1 = ct.ions[ion2].iyend - ct.ions[ion2].iystart + 1;
    izz1 = ct.ions[ion2].izend - ct.ions[ion2].izstart + 1;
    incx1 = iyy1 * izz1;
    incy1 = izz1;

    for (ix = xlow1; ix <= xhigh1; ix++)
    {
        ix1 = (ix - state1->ixmin) * incx;
        ix2 = (ix - ct.ions[ion2].ixstart) * incx1;

        for (iy = ylow1; iy <= yhigh1; iy++)
        {
            iy1 = (iy - state1->iymin) * incy;
            iy2 = (iy - ct.ions[ion2].iystart) * incy1;

            for (iz = zlow1; iz <= zhigh1; iz++)
            {
                iz1 = iz - state1->izmin;
                iz2 = iz - ct.ions[ion2].izstart;

                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;

                work[idx1] += pptr[idx2];
            }                   /* end for iz = zlow1 */

            for (iz = zlow2; iz <= zhigh2; iz++)
            {
                iz1 = iz - state1->izmin;
                iz2 = iz - ct.ions[ion2].izstart - zshift;

                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;

                work[idx1] += pptr[idx2];
            }                   /* end for iz = zlow2 */
        }                       /* end for iy = ylow1 */

        for (iy = ylow2; iy <= yhigh2; iy++)
        {
            iy1 = (iy - state1->iymin) * incy;
            iy2 = (iy - ct.ions[ion2].iystart - yshift) * incy1;

            for (iz = zlow1; iz <= zhigh1; iz++)
            {
                iz1 = iz - state1->izmin;
                iz2 = iz - ct.ions[ion2].izstart;

                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;

                work[idx1] += pptr[idx2];
            }                   /* end for iz = zlow1 */

            for (iz = zlow2; iz <= zhigh2; iz++)
            {
                iz1 = iz - state1->izmin;
                iz2 = iz - ct.ions[ion2].izstart - zshift;

                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;

                work[idx1] += pptr[idx2];
            }                   /* end for iz = zlow2 */
        }                       /* end for iy = ylow2 */
    }                           /* end for ix = xlow1 */

    for (ix = xlow2; ix <= xhigh2; ix++)
    {
        ix1 = (ix - state1->ixmin) * incx;
        ix2 = (ix - ct.ions[ion2].ixstart - xshift) * incx1;

        for (iy = ylow1; iy <= yhigh1; iy++)
        {
            iy1 = (iy - state1->iymin) * incy;
            iy2 = (iy - ct.ions[ion2].iystart) * incy1;

            for (iz = zlow1; iz <= zhigh1; iz++)
            {
                iz1 = iz - state1->izmin;
                iz2 = iz - ct.ions[ion2].izstart;

                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;

                work[idx1] += pptr[idx2];
            }                   /* end for iz = zlow1 */

            for (iz = zlow2; iz <= zhigh2; iz++)
            {
                iz1 = iz - state1->izmin;
                iz2 = iz - ct.ions[ion2].izstart - zshift;

                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;

                work[idx1] += pptr[idx2];
            }                   /* end for iz = zlow2 */
        }                       /* end for iy = ylow1 */

        for (iy = ylow2; iy <= yhigh2; iy++)
        {
            iy1 = (iy - state1->iymin) * incy;
            iy2 = (iy - ct.ions[ion2].iystart - yshift) * incy1;

            for (iz = zlow1; iz <= zhigh1; iz++)
            {
                iz1 = iz - state1->izmin;
                iz2 = iz - ct.ions[ion2].izstart;

                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;

                work[idx1] += pptr[idx2];
            }                   /* end for iz = zlow1 */

            for (iz = zlow2; iz <= zhigh2; iz++)
            {
                iz1 = iz - state1->izmin;
                iz2 = iz - ct.ions[ion2].izstart - zshift;

                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;

                work[idx1] += pptr[idx2];
            }                   /* end for iz = zlow2 */
        }                       /* end for iy = ylow2 */
    }                           /* end for ix = xlow2 */


}
