/************************** SVN Revision Information **************************
 **    $Id: theta_phi_new.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
/*
	calculate res_i = theta_ij * phi_j
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void theta_phi_new(int st1, int st2, REAL theta_ion, REAL * st2_psi,
                   REAL * state1_psi, int mode, STATE * states)
{
    int xlow1, xhigh1, xlow2, xhigh2, xshift;
    int ylow1, yhigh1, ylow2, yhigh2, yshift;
    int zlow1, zhigh1, zlow2, zhigh2, zshift;
    int iyy, izz, iyy1, izz1;
    int incx, incy, incx1, incy1;
    int ix, iy, iz, ix1, ix2, iy1, iy2, iz1, iz2;
    int idx1, idx2;
    int index, xshift1, xshift2, yshift1, yshift2, zshift1, zshift2;
    int zlength1, zlength2;

    if (mode == 0)
        index = (st1 - ct.state_begin) * ct.num_states + st2;
    if (mode == 1)
        index = (st2 - ct.state_begin) * ct.num_states + st1;
    if (mode < 0 || mode > 1)
        error_handler(" mode is not right");


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

    iyy = states[st1].iymax - states[st1].iymin + 1;
    izz = states[st1].izmax - states[st1].izmin + 1;
    incx = iyy * izz;
    incy = izz;

    iyy1 = states[st2].iymax - states[st2].iymin + 1;
    izz1 = states[st2].izmax - states[st2].izmin + 1;
    incx1 = iyy1 * izz1;
    incy1 = izz1;

    for (ix = xlow1; ix <= xhigh1; ix++)
    {
        ix1 = (ix - states[st1].ixmin) * incx;
        ix2 = (ix - states[st2].ixmin) * incx1;

        for (iy = ylow1; iy <= yhigh1; iy++)
        {
            iy1 = (iy - states[st1].iymin) * incy;
            iy2 = (iy - states[st2].iymin) * incy1;
            idx1 = ix1 + iy1 + zlow1 - states[st1].izmin;
            idx2 = ix2 + iy2 + zlow1 - states[st2].izmin;
            for (iz = 0; iz < zlength1; iz++)
                state1_psi[idx1 + iz] += theta_ion * st2_psi[idx2 + iz];

            idx1 = ix1 + iy1 + zlow2 - states[st1].izmin - zshift1;
            idx2 = ix2 + iy2 + zlow2 - states[st2].izmin - zshift2;
            for (iz = 0; iz < zlength2; iz++)
                state1_psi[idx1 + iz] += theta_ion * st2_psi[idx2 + iz];
        }

        for (iy = ylow2; iy <= yhigh2; iy++)
        {
            iy1 = (iy - states[st1].iymin - yshift1) * incy;
            iy2 = (iy - states[st2].iymin - yshift2) * incy1;
            idx1 = ix1 + iy1 + zlow1 - states[st1].izmin;
            idx2 = ix2 + iy2 + zlow1 - states[st2].izmin;
            for (iz = 0; iz < zlength1; iz++)
                state1_psi[idx1 + iz] += theta_ion * st2_psi[idx2 + iz];

            idx1 = ix1 + iy1 + zlow2 - states[st1].izmin - zshift1;
            idx2 = ix2 + iy2 + zlow2 - states[st2].izmin - zshift2;
            for (iz = 0; iz < zlength2; iz++)
                state1_psi[idx1 + iz] += theta_ion * st2_psi[idx2 + iz];
        }
    }

    for (ix = xlow2; ix <= xhigh2; ix++)
    {
        ix1 = (ix - states[st1].ixmin - xshift1) * incx;
        ix2 = (ix - states[st2].ixmin - xshift2) * incx1;
        for (iy = ylow1; iy <= yhigh1; iy++)
        {
            iy1 = (iy - states[st1].iymin) * incy;
            iy2 = (iy - states[st2].iymin) * incy1;
            idx1 = ix1 + iy1 + zlow1 - states[st1].izmin;
            idx2 = ix2 + iy2 + zlow1 - states[st2].izmin;
            for (iz = 0; iz < zlength1; iz++)
                state1_psi[idx1 + iz] += theta_ion * st2_psi[idx2 + iz];

            idx1 = ix1 + iy1 + zlow2 - states[st1].izmin - zshift1;
            idx2 = ix2 + iy2 + zlow2 - states[st2].izmin - zshift2;
            for (iz = 0; iz < zlength2; iz++)
                state1_psi[idx1 + iz] += theta_ion * st2_psi[idx2 + iz];
        }

        for (iy = ylow2; iy <= yhigh2; iy++)
        {
            iy1 = (iy - states[st1].iymin - yshift1) * incy;
            iy2 = (iy - states[st2].iymin - yshift2) * incy1;
            idx1 = ix1 + iy1 + zlow1 - states[st1].izmin;
            idx2 = ix2 + iy2 + zlow1 - states[st2].izmin;
            for (iz = 0; iz < zlength1; iz++)
                state1_psi[idx1 + iz] += theta_ion * st2_psi[idx2 + iz];

            idx1 = ix1 + iy1 + zlow2 - states[st1].izmin - zshift1;
            idx2 = ix2 + iy2 + zlow2 - states[st2].izmin - zshift2;
            for (iz = 0; iz < zlength2; iz++)
                state1_psi[idx1 + iz] += theta_ion * st2_psi[idx2 + iz];
        }
    }
}
