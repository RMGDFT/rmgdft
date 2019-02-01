/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*


dot_product_orbit_nl.c

dot_product of (orbit,  non-local projector )

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "blas.h"
#include "Kbpsi.h"

static void inline inline_betapsi(int zlength, int num_proj, double *prj, 
        double *psi, double *kbpsi)
{
    for(int i = 0; i < num_proj; i++)
        for(int j = 0; j < zlength; j++)
            kbpsi[i] += get_vel() * prj[i * ct.max_nlpoints + j] * psi[j];

    //double alpha = get_vel();
    //double one = 1.0;
    //int ione = 1;
    //dgemv("T", &zlength, &num_proj, &alpha, prj, &ct.max_nlpoints, 
    //        psi, &ione, &one, kbpsi, &ione); 
}
void DotProductOrbitNl(STATE *st1, int ion2, double *
        psi, double * prjptr, ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl,
        int num_proj, double *kbpsi)
{

    int xlow1, xhigh1, xlow2, xhigh2, xshift;
    int ylow1, yhigh1, ylow2, yhigh2, yshift;
    int zlow1, zhigh1, zlow2, zhigh2, zshift;
    int iyy, izz, iyy1, izz1;
    int incx, incy, incx1, incy1;
    int ix, iy, ix1, ix2, iy1, iy2, iz1, iz2;
    int idx1, idx2;
    int index;
    int zlength1, zlength2;

    index = (st1->index - ct.state_begin) * ct.num_ions + ion2;

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

    zlength1 = zhigh1 - zlow1 +1;
    zlength2 = zhigh2 - zlow2 +1;

    iyy = st1->iymax - st1->iymin + 1;
    izz = st1->izmax - st1->izmin + 1;
    incx = iyy * izz;
    incy = izz;

    iyy1 = ct.ions[ion2].iyend - ct.ions[ion2].iystart + 1;
    izz1 = ct.ions[ion2].izend - ct.ions[ion2].izstart + 1;
    incx1 = iyy1 * izz1;
    incy1 = izz1;

    for(int i = 0; i < num_proj; i++) kbpsi[i] = 0.0;

    if(zlength1 > 0)
    {
        iz1 = zlow1 - st1->izmin;
        iz2 = zlow1 - ct.ions[ion2].izstart;
        for (ix = xlow1; ix <= xhigh1; ix++)
        {
            ix1 = (ix - st1->ixmin) * incx;
            ix2 = (ix - ct.ions[ion2].ixstart) * incx1;

            for (iy = ylow1; iy <= yhigh1; iy++)
            {
                iy1 = (iy - st1->iymin) * incy;
                iy2 = (iy - ct.ions[ion2].iystart) * incy1;

                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                inline_betapsi(zlength1, num_proj, &prjptr[idx2], &psi[idx1], kbpsi); 
            }


            for (iy = ylow2; iy <= yhigh2; iy++)
            {
                iy1 = (iy - st1->iymin) * incy;
                iy2 = (iy - ct.ions[ion2].iystart - yshift) * incy1;

                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                inline_betapsi(zlength1, num_proj, &prjptr[idx2], &psi[idx1], kbpsi); 
            }
        }                           /* end for ix = xlow1 */

        for (ix = xlow2; ix <= xhigh2; ix++)
        {
            ix1 = (ix - st1->ixmin) * incx;
            ix2 = (ix - ct.ions[ion2].ixstart - xshift) * incx1;

            for (iy = ylow1; iy <= yhigh1; iy++)
            {
                iy1 = (iy - st1->iymin) * incy;
                iy2 = (iy - ct.ions[ion2].iystart) * incy1;
                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                inline_betapsi(zlength1, num_proj, &prjptr[idx2], &psi[idx1], kbpsi); 
            }                       /* end for iy = ylow1 */

            for (iy = ylow2; iy <= yhigh2; iy++)
            {
                iy1 = (iy - st1->iymin) * incy;
                iy2 = (iy - ct.ions[ion2].iystart - yshift) * incy1;
                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                inline_betapsi(zlength1, num_proj, &prjptr[idx2], &psi[idx1], kbpsi); 
            }
        }                           /* end for ix = xlow2 */


    }

    if(zlength2 > 0)
    {
        iz1 = zlow2 - st1->izmin;
        iz2 = zlow2 - ct.ions[ion2].izstart-zshift;
        for (ix = xlow1; ix <= xhigh1; ix++)
        {
            ix1 = (ix - st1->ixmin) * incx;
            ix2 = (ix - ct.ions[ion2].ixstart) * incx1;

            for (iy = ylow1; iy <= yhigh1; iy++)
            {
                iy1 = (iy - st1->iymin) * incy;
                iy2 = (iy - ct.ions[ion2].iystart) * incy1;

                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                inline_betapsi(zlength2, num_proj, &prjptr[idx2], &psi[idx1], kbpsi); 
            }

            for (iy = ylow2; iy <= yhigh2; iy++)
            {
                iy1 = (iy - st1->iymin) * incy;
                iy2 = (iy - ct.ions[ion2].iystart - yshift) * incy1;

                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                inline_betapsi(zlength2, num_proj, &prjptr[idx2], &psi[idx1], kbpsi); 
            }
        }                           /* end for ix = xlow1 */

        for (ix = xlow2; ix <= xhigh2; ix++)
        {
            ix1 = (ix - st1->ixmin) * incx;
            ix2 = (ix - ct.ions[ion2].ixstart - xshift) * incx1;

            for (iy = ylow1; iy <= yhigh1; iy++)
            {
                iy1 = (iy - st1->iymin) * incy;
                iy2 = (iy - ct.ions[ion2].iystart) * incy1;
                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                inline_betapsi(zlength2, num_proj, &prjptr[idx2], &psi[idx1], kbpsi); 
            }                       /* end for iy = ylow1 */

            for (iy = ylow2; iy <= yhigh2; iy++)
            {
                iy1 = (iy - st1->iymin) * incy;
                iy2 = (iy - ct.ions[ion2].iystart - yshift) * incy1;
                idx1 = ix1 + iy1 + iz1;
                idx2 = ix2 + iy2 + iz2;
                inline_betapsi(zlength2, num_proj, &prjptr[idx2], &psi[idx1], kbpsi); 
            }
        }                           /* end for ix = xlow2 */

    }

}
