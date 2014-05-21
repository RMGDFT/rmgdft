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
#include "prototypes_on.h"
#include "init_var.h"



void dot_product_orbit_xyz_orbit(STATE *orbit1, STATE *orbit2, double *X0, double *Y0, double *Z0)
{

    int xlow1, xhigh1, xlow2, xhigh2, xshift;
    int ylow1, yhigh1, ylow2, yhigh2, yshift;
    int zlow1, zhigh1, zlow2, zhigh2, zshift;
    int iyy, izz, iyy1, izz1;
    int incx, incy, incx1, incy1;
    int zcase, ix, iy, iz;
    int ix1, ix2, iy1, iy2, idx1, idx2;
    double time1;
    int index, xshift1, xshift2, yshift1, yshift2, zshift1, zshift2;
    int zlength1, zlength2;
    double *psi1, *psi2;
    double hxgrid, hygrid, hzgrid;

    int ione = 1;


    hxgrid = get_hxgrid() * get_xside();
    hygrid = get_hxgrid() * get_yside();
    hzgrid = get_hxgrid() * get_zside();

    time1 = my_crtc();

    *X0 = 0.0;
    *Y0 = 0.0;
    *Z0 = 0.0;

    if (orbit1->index < ct.state_begin || orbit1->index >= ct.state_end)
        error_handler("orbit1 is not in this PE");
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
    psi1 = orbit1->psiR;
    psi2 = orbit2->psiR;

    if(zlength1 >0 && zlength2 >0) 
    {
        zcase = 1;
    }
    else if( zlength1 >0)
    {
        zcase = 2;
    }
    else if(zlength2 >0)
    {
        zcase = 3;
    }
    else
    {
        printf("\n zlength1= %d  zlength2=%d\n", zlength1, zlength2);
        printf("orbit1 %d  orbit2 %d has no overlap\n", orbit1->index, orbit2->index);
        exit(0);
    }


    switch(zcase)
    {

        case 1:
            break;
            for (ix = xlow1; ix <= xhigh1; ix++)
            {
                ix1 = (ix - orbit1->ixmin) * incx;
                ix2 = (ix - orbit2->ixmin) * incx1;

                for (iy = ylow1; iy <= yhigh1; iy++)
                {
                    iy1 = (iy - orbit1->iymin) * incy;
                    iy2 = (iy - orbit2->iymin) * incy1;

                    for (iz = zlow1; iz <= zhigh1; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

                    for (iz = zlow2; iz <= zhigh2; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin - zshift1;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin - zshift2;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }


                }

                for (iy = ylow2; iy <= yhigh2; iy++)
                {
                    iy1 = (iy - orbit1->iymin - yshift1) * incy;
                    iy2 = (iy - orbit2->iymin - yshift2) * incy1;

                    for (iz = zlow1; iz <= zhigh1; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

                    for (iz = zlow2; iz <= zhigh2; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin - zshift1;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin - zshift2;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }


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

                    for (iz = zlow1; iz <= zhigh1; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

                    for (iz = zlow2; iz <= zhigh2; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin - zshift1;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin - zshift2;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }


                }

                for (iy = ylow2; iy <= yhigh2; iy++)
                {
                    iy1 = (iy - orbit1->iymin - yshift1) * incy;
                    iy2 = (iy - orbit2->iymin - yshift2) * incy1;

                    for (iz = zlow1; iz <= zhigh1; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

                    for (iz = zlow2; iz <= zhigh2; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin - zshift1;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin - zshift2;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

                }
            }

            break;

        case 2:
            for (ix = xlow1; ix <= xhigh1; ix++)
            {
                ix1 = (ix - orbit1->ixmin) * incx;
                ix2 = (ix - orbit2->ixmin) * incx1;

                for (iy = ylow1; iy <= yhigh1; iy++)
                {
                    iy1 = (iy - orbit1->iymin) * incy;
                    iy2 = (iy - orbit2->iymin) * incy1;

                    for (iz = zlow1; iz <= zhigh1; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

                }

                for (iy = ylow2; iy <= yhigh2; iy++)
                {
                    iy1 = (iy - orbit1->iymin - yshift1) * incy;
                    iy2 = (iy - orbit2->iymin - yshift2) * incy1;

                    for (iz = zlow1; iz <= zhigh1; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

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

                    for (iz = zlow1; iz <= zhigh1; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

                }

                for (iy = ylow2; iy <= yhigh2; iy++)
                {
                    iy1 = (iy - orbit1->iymin - yshift1) * incy;
                    iy2 = (iy - orbit2->iymin - yshift2) * incy1;

                    for (iz = zlow1; iz <= zhigh1; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

                }
            }


            break;

        case 3:
            break;
            for (ix = xlow1; ix <= xhigh1; ix++)
            {
                ix1 = (ix - orbit1->ixmin) * incx;
                ix2 = (ix - orbit2->ixmin) * incx1;

                for (iy = ylow1; iy <= yhigh1; iy++)
                {
                    iy1 = (iy - orbit1->iymin) * incy;
                    iy2 = (iy - orbit2->iymin) * incy1;

                    for (iz = zlow2; iz <= zhigh2; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin -zshift1;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin -zshift2;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

                }

                for (iy = ylow2; iy <= yhigh2; iy++)
                {
                    iy1 = (iy - orbit1->iymin - yshift1) * incy;
                    iy2 = (iy - orbit2->iymin - yshift2) * incy1;

                    for (iz = zlow2; iz <= zhigh2; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin -zshift1;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin -zshift2;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

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

                    for (iz = zlow2; iz <= zhigh2; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

                }

                for (iy = ylow2; iy <= yhigh2; iy++)
                {
                    iy1 = (iy - orbit1->iymin - yshift1) * incy;
                    iy2 = (iy - orbit2->iymin - yshift2) * incy1;

                    for (iz = zlow2; iz <= zhigh2; iz++)
                    {

                        idx1 = ix1 + iy1 + iz - orbit1->izmin -zshift1;
                        idx2 = ix2 + iy2 + iz - orbit2->izmin -zshift2;
                        *X0 += psi1[idx1] * psi2[idx2] * ix * hxgrid;
                        *Y0 += psi1[idx1] * psi2[idx2] * iy * hygrid;
                        *Z0 += psi1[idx1] * psi2[idx2] * iz * hzgrid;
                    }

                }
            }

            break;

    }

    time1 = my_crtc() - time1;
    rmg_timings(DOT_PRODUCT, time1);
}

