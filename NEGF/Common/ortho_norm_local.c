/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var_negf.h"
#include "LCR.h"

#define 	MAX_ORBIT_ON_ION 	24

void state_minus_state (STATE *state1, STATE *state2, rmg_double_t factor);


void ortho_norm_local (STATE *states)
{
    int ione = 1, i, j, k;
    int n;
    double time1, time2;
    int num_state_on_this_ion;
    rmg_double_t norm;
    rmg_double_t overlap[MAX_ORBIT_ON_ION];

    time1 = my_crtc ();

    if (pct.gridpe == 0)
        printf ("\n LOCAL ORTHONORMALIZATION  ");

    for (k = ct.ion_begin; k < ct.ion_end; k++)
    {
        num_state_on_this_ion = 0;
        for (j = ct.state_begin; j < ct.state_end; j++)
        {
            if (state_to_ion[j] == k)
            {
                for (n = 0; n < MAX_ORBIT_ON_ION; n++)
                    overlap[n] = 0.0;
                num_state_on_this_ion++;
                for (i = 1; i < num_state_on_this_ion; i++)
                {
                    assert (states[j].radius >= states[j - i].radius);
                    overlap[i] = dot_product_orbit_orbit (&states[j], &states[j - i]);
                }
                for (i = 1; i < num_state_on_this_ion; i++)
                {
                    state_minus_state (&states[j], &states[j - i], overlap[i]);
                }
                norm = QMD_ddot (states[j].size, states[j].psiR, ione, states[j].psiR, ione);
                norm = 1.0 / (sqrt (norm));
                sscal (&states[j].size, &norm, states[j].psiR, &ione);

            }
        }
    }

    for (j = ct.state_begin; j < ct.state_end; j++)
    {
        norm = 1.0 / (sqrt (ct.vel));
        sscal (&states[j].size, &norm, states[j].psiR, &ione);
    }

    time2 = my_crtc ();
    rmg_timings (ORTHONORM_TIME, (time2 - time1));
}


void state_minus_state (STATE *orbit1, STATE *orbit2, rmg_double_t factor)
{
    int xlow1, xhigh1, xlow2, xhigh2, xshift;
    int ylow1, yhigh1, ylow2, yhigh2, yshift;
    int zlow1, zhigh1, zlow2, zhigh2, zshift;
    int iyy, izz, iyy1, izz1;
    int incx, incy, incx1, incy1;
    int ix, iy, iz;
    int ix1, ix2, iy1, iy2, idx1, idx2;
    int index, xshift1, xshift2, yshift1, yshift2, zshift1, zshift2;
    int zlength1, zlength2;
    rmg_double_t *p1, *p2;
    rmg_double_t *psi1, *psi2;
    int mode;

    mode = 0;
    index = (orbit1->index - ct.state_begin) * ct.num_states + orbit2->index;

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

    xshift1 = mode * xshift;
    xshift2 = (1 - mode) * xshift;
    yshift1 = mode * yshift;
    yshift2 = (1 - mode) * yshift;
    zshift1 = mode * zshift;
    zshift2 = (1 - mode) * zshift;

    iyy = orbit1->iymax - orbit1->iymin + 1;
    izz = orbit1->izmax - orbit1->izmin + 1;
    incx = iyy * izz;
    incy = izz;

    iyy1 = orbit2->iymax - orbit2->iymin + 1;
    izz1 = orbit2->izmax - orbit2->izmin + 1;
    incx1 = iyy1 * izz1;
    incy1 = izz1;
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
                    p1[iz] = p1[iz] - factor * p2[iz];
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
                    p1[iz] = p1[iz] - factor * p2[iz];
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
                    p1[iz] = p1[iz] - factor * p2[iz];
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
                    p1[iz] = p1[iz] - factor * p2[iz];
            }
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
                    p1[iz] = p1[iz] - factor * p2[iz];
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
                    p1[iz] = p1[iz] - factor * p2[iz];
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
                idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                p1 = &psi1[idx1];
                p2 = &psi2[idx2];
                for (iz = 0; iz < zlength2; iz++)
                    p1[iz] = p1[iz] - factor * p2[iz];
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
                    p1[iz] = p1[iz] - factor * p2[iz];
            }
        }
    }                           /* end zlength2 */

}
