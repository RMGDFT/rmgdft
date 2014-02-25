/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
	dot_product of (orbit, orbit)

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "init_var.h"

#define   SUM_SIZE 	1000

rmg_double_t dot_product_orbit_orbit (STATE *orbit1, STATE *orbit2)
{

    int xlow1, xhigh1, xlow2, xhigh2, xshift;
    int ylow1, yhigh1, ylow2, yhigh2, yshift;
    int zlow1, zhigh1, zlow2, zhigh2, zshift;
    int iyy, izz, iyy1, izz1;
    int incx, incy, incx1, incy1;
    int ix, iy, iz;
    int ix1, ix2, iy1, iy2, idx1, idx2;
    rmg_double_t time1;
    int index, xshift1, xshift2, yshift1, yshift2, zshift1, zshift2;
    int zlength1, zlength2;
    rmg_double_t z_sum[SUM_SIZE];
    rmg_double_t *p1, *p2;
    rmg_double_t *psi1, *psi2;
    rmg_double_t result;


    time1 = my_crtc ();

    assert (max (orbit1->orbit_nz, orbit2->orbit_nz) < SUM_SIZE);
    if (orbit1->index < ct.state_begin || orbit1->index >= ct.state_end)
    {
        printf ("\n orbit1-> %d, orbit2-> %d ", orbit1->index, orbit2->index);
        fflush (NULL);
        error_handler ("orbit1->is not in this PE");
    }
    index = (orbit1->index - ct.state_begin) * ct.num_states + orbit2->index;

    result = 0.0;
    xlow1 = orbit_overlap_region[index].xlow1;
    xhigh1 = orbit_overlap_region[index].xhigh1;
    xlow2 = orbit_overlap_region[index].xlow2;
    xhigh2 = orbit_overlap_region[index].xhigh2;
    xshift = orbit_overlap_region[index].xshift;

    ylow1 = orbit_overlap_region[index].ylow1;
    yhigh1 = orbit_overlap_region[index].yhigh1;
    ylow2 = orbit_overlap_region[index].ylow2;
    yhigh2 = orbit_overlap_region[index].yhigh2;
    yshift = orbit_overlap_region[index].yshift;

    zlow1 = orbit_overlap_region[index].zlow1;
    zhigh1 = orbit_overlap_region[index].zhigh1;
    zlow2 = orbit_overlap_region[index].zlow2;
    zhigh2 = orbit_overlap_region[index].zhigh2;
    zshift = orbit_overlap_region[index].zshift;

    zlength1 = zhigh1 - zlow1 + 1;
    zlength2 = zhigh2 - zlow2 + 1;

    xshift1 = 0.0;
    xshift2 = xshift;
    yshift1 = 0.0;
    yshift2 = yshift;
    zshift1 = 0.0;
    zshift2 = zshift;

    iyy = orbit1->iymax - orbit1->iymin + 1;
    izz = orbit1->izmax - orbit1->izmin + 1;
    incx = iyy * izz;
    incy = izz;

    iyy1 = orbit2->iymax - orbit2->iymin + 1;
    izz1 = orbit2->izmax - orbit2->izmin + 1;
    incx1 = iyy1 * izz1;
    incy1 = izz1;
    for (iz = 0; iz < SUM_SIZE; iz++)
        z_sum[iz] = 0.0;
    psi1 = orbit1->psiR;
    psi2 = orbit2->psiR;

    if (zlength1 > 0)
    {
        for (ix = xlow1; ix <= xhigh1; ix++)
        {
            ix1 = (ix - orbit1->ixmin) * incx;
            ix2 = (ix - orbit2->ixmin) * incx1;

            for (iy = ylow1; iy <= yhigh1; iy++)
            {
                iy1 = (iy - orbit1->iymin) * incy;
                iy2 = (iy - orbit2->iymin) * incy1;
                idx1 = ix1 + iy1 + zlow1 - orbit1->izmin;
                idx2 = ix2 + iy2 + zlow1 - orbit2->izmin;
                p1 = &psi1[idx1];
                p2 = &psi2[idx2];
                for (iz = 0; iz < zlength1; iz++)
                    z_sum[iz] += p1[iz] * p2[iz];
            }

            for (iy = ylow2; iy <= yhigh2; iy++)
            {
                iy1 = (iy - orbit1->iymin - yshift1) * incy;
                iy2 = (iy - orbit2->iymin - yshift2) * incy1;
                idx1 = ix1 + iy1 + zlow1 - orbit1->izmin;
                idx2 = ix2 + iy2 + zlow1 - orbit2->izmin;
                p1 = &psi1[idx1];
                p2 = &psi2[idx2];
                for (iz = 0; iz < zlength1; iz++)
                    z_sum[iz] += p1[iz] * p2[iz];
            }
        }

        for (ix = xlow2; ix <= xhigh2; ix++)
        {
            ix1 = (ix - orbit1->ixmin - xshift1) * incx;
            ix2 = (ix - orbit2->ixmin - xshift2) * incx1;
            for (iy = ylow1; iy <= yhigh1; iy++)
            {
                iy1 = (iy - orbit1->iymin) * incy;
                iy2 = (iy - orbit2->iymin) * incy1;
                idx1 = ix1 + iy1 + zlow1 - orbit1->izmin;
                idx2 = ix2 + iy2 + zlow1 - orbit2->izmin;
                p1 = &psi1[idx1];
                p2 = &psi2[idx2];
                for (iz = 0; iz < zlength1; iz++)
                    z_sum[iz] += p1[iz] * p2[iz];
            }

            for (iy = ylow2; iy <= yhigh2; iy++)
            {
                iy1 = (iy - orbit1->iymin - yshift1) * incy;
                iy2 = (iy - orbit2->iymin - yshift2) * incy1;
                idx1 = ix1 + iy1 + zlow1 - orbit1->izmin;
                idx2 = ix2 + iy2 + zlow1 - orbit2->izmin;
                p1 = &psi1[idx1];
                p2 = &psi2[idx2];
                for (iz = 0; iz < zlength1; iz++)
                    z_sum[iz] += p1[iz] * p2[iz];
            }
        }
        for (iz = 0; iz < zlength1; iz++)
        {
            result += z_sum[iz];
            z_sum[iz] = 0.0;
        }
    }                           /* end zlength1 */




    if (zlength2 > 0)
    {
        for (ix = xlow1; ix <= xhigh1; ix++)
        {
            ix1 = (ix - orbit1->ixmin) * incx;
            ix2 = (ix - orbit2->ixmin) * incx1;
            for (iy = ylow1; iy <= yhigh1; iy++)
            {
                iy1 = (iy - orbit1->iymin) * incy;
                iy2 = (iy - orbit2->iymin) * incy1;
                idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                p1 = &psi1[idx1];
                p2 = &psi2[idx2];
                for (iz = 0; iz < zlength2; iz++)
                    z_sum[iz] += p1[iz] * p2[iz];
            }

            for (iy = ylow2; iy <= yhigh2; iy++)
            {
                iy1 = (iy - orbit1->iymin - yshift1) * incy;
                iy2 = (iy - orbit2->iymin - yshift2) * incy1;
                idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                p1 = &psi1[idx1];
                p2 = &psi2[idx2];
                for (iz = 0; iz < zlength2; iz++)
                    z_sum[iz] += p1[iz] * p2[iz];
            }
        }

        for (ix = xlow2; ix <= xhigh2; ix++)
        {
            ix1 = (ix - orbit1->ixmin - xshift1) * incx;
            ix2 = (ix - orbit2->ixmin - xshift2) * incx1;
            for (iy = ylow1; iy <= yhigh1; iy++)
            {
                iy1 = (iy - orbit1->iymin) * incy;
                iy2 = (iy - orbit2->iymin) * incy1;
                idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                idx2 = ix2 + iy2 + zlow2 - orbit1->izmin - zshift2;
                p1 = &psi1[idx1];
                p2 = &psi2[idx2];
                for (iz = 0; iz < zlength2; iz++)
                    z_sum[iz] += p1[iz] * p2[iz];
            }

            for (iy = ylow2; iy <= yhigh2; iy++)
            {
                iy1 = (iy - orbit1->iymin - yshift1) * incy;
                iy2 = (iy - orbit2->iymin - yshift2) * incy1;
                idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                p1 = &psi1[idx1];
                p2 = &psi2[idx2];
                for (iz = 0; iz < zlength2; iz++)
                    z_sum[iz] += p1[iz] * p2[iz];
            }
        }

        for (iz = 0; iz < zlength2; iz++)
        {
            result += z_sum[iz];
            z_sum[iz] = 0.0;
        }

    }                           /* end zlength2 */

    time1 = my_crtc () - time1;
    rmg_timings (DOT_PRODUCT, time1);
    return result;
}
