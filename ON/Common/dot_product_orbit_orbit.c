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


void dot_product_orbit_orbit(STATE *orbit1, STATE *orbit2, STATE
*orbit3, double *H, double *S)
{

    int xlow1, xhigh1, xlow2, xhigh2, xshift;
    int ylow1, yhigh1, ylow2, yhigh2, yshift;
    int zlow1, zhigh1, zlow2, zhigh2, zshift;
    int iyy, izz, iyy1, izz1;
    int incx, incy, incx1, incy1;
    int zcase, ix, iy, iz;
    int ix1, ix2, iy1, iy2, idx1, idx2;
    REAL time1;
    int index, xshift1, xshift2, yshift1, yshift2, zshift1, zshift2;
    int zlength1, zlength2;
    REAL *p1, *p2, *p3;
    REAL *psi1, *psi2, *psi3;

    int ione = 1;


    time1 = my_crtc();

    *H = 0.0;
    *S = 0.0;

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
    psi3 = orbit3->psiR;

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
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += ddot(&zlength1, p3, &ione, p2, &ione);

                    idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                    idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += ddot(&zlength2, p3, &ione, p2, &ione);

                }

                for (iy = ylow2; iy <= yhigh2; iy++)
                {
                    iy1 = (iy - orbit1->iymin - yshift1) * incy;
                    iy2 = (iy - orbit2->iymin - yshift2) * incy1;

                    idx1 = ix1 + iy1 + zlow1 - orbit1->izmin;
                    idx2 = ix2 + iy2 + zlow1 - orbit2->izmin;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += ddot(&zlength1, p3, &ione, p2, &ione);

                    idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                    idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += ddot(&zlength2, p3, &ione, p2, &ione);

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
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += ddot(&zlength1, p3, &ione, p2, &ione);

                    idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                    idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += ddot(&zlength2, p3, &ione, p2, &ione);

                }

                for (iy = ylow2; iy <= yhigh2; iy++)
                {
                    iy1 = (iy - orbit1->iymin - yshift1) * incy;
                    iy2 = (iy - orbit2->iymin - yshift2) * incy1;
                    idx1 = ix1 + iy1 + zlow1 - orbit1->izmin;
                    idx2 = ix2 + iy2 + zlow1 - orbit2->izmin;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += ddot(&zlength1, p3, &ione, p2, &ione);

                    idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                    idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += ddot(&zlength2, p3, &ione, p2, &ione);

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
                    idx1 = ix1 + iy1 + zlow1 - orbit1->izmin;
                    idx2 = ix2 + iy2 + zlow1 - orbit2->izmin;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += ddot(&zlength1, p3, &ione, p2, &ione);

                }

                for (iy = ylow2; iy <= yhigh2; iy++)
                {
                    iy1 = (iy - orbit1->iymin - yshift1) * incy;
                    iy2 = (iy - orbit2->iymin - yshift2) * incy1;

                    idx1 = ix1 + iy1 + zlow1 - orbit1->izmin;
                    idx2 = ix2 + iy2 + zlow1 - orbit2->izmin;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += ddot(&zlength1, p3, &ione, p2, &ione);

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
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += ddot(&zlength1, p3, &ione, p2, &ione);

                }

                for (iy = ylow2; iy <= yhigh2; iy++)
                {
                    iy1 = (iy - orbit1->iymin - yshift1) * incy;
                    iy2 = (iy - orbit2->iymin - yshift2) * incy1;
                    idx1 = ix1 + iy1 + zlow1 - orbit1->izmin;
                    idx2 = ix2 + iy2 + zlow1 - orbit2->izmin;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += ddot(&zlength1, p3, &ione, p2, &ione);

                }
            }


            break;

        case 3:
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
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += ddot(&zlength2, p3, &ione, p2, &ione);

                }

                for (iy = ylow2; iy <= yhigh2; iy++)
                {
                    iy1 = (iy - orbit1->iymin - yshift1) * incy;
                    iy2 = (iy - orbit2->iymin - yshift2) * incy1;
                    idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                    idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += ddot(&zlength2, p3, &ione, p2, &ione);

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
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += ddot(&zlength2, p3, &ione, p2, &ione);

                }

                for (iy = ylow2; iy <= yhigh2; iy++)
                {
                    iy1 = (iy - orbit1->iymin - yshift1) * incy;
                    iy2 = (iy - orbit2->iymin - yshift2) * incy1;
                    idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                    idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += ddot(&zlength2, p3, &ione, p2, &ione);

                }
            }

            break;

    }

    time1 = my_crtc() - time1;
    rmg_timings(DOT_PRODUCT, time1);
}

