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

static double inline inline_ddot(int *length, double *p1, int *ione, double *p2, int *iione)
{
    double sum = 0.0;
    for(int idx=0;idx < *length;idx++) sum += p1[idx] * p2[idx];
    return sum;
}


void DotProductOrbitOrbit(STATE *orbit1, STATE *orbit2, STATE
        *orbit3, double *H, double *S, ORBITAL_PAIR *onepair)
{

    int xlow1, xhigh1, xlow2, xhigh2, xshift;
    int ylow1, yhigh1, ylow2, yhigh2, yshift;
    int zlow1, zhigh1, zlow2, zhigh2, zshift;
    int iyy, izz, iyy1, izz1;
    int incx, incy, incx1, incy1;
    int zcase, ix, iy;
    int ix1, ix2, iy1, iy2, idx1, idx2;
    int index, xshift1, xshift2, yshift1, yshift2, zshift1, zshift2;
    int zlength1, zlength2;
    double *p1, *p2, *p3;
    double *psi1, *psi2, *psi3;

    int ione = 1;


    void *RT1 = NULL;
    if(ct.verbose) RT1 = new RmgTimer("4-get_HS: orbit_dot_orbit: dotproduct");

    *H = 0.0;
    *S = 0.0;

    if (orbit1->index < ct.state_begin || orbit1->index >= ct.state_end)
        error_handler("orbit1 is not in this PE");

    xlow1 = onepair->xlow1;
    xhigh1 = onepair->xhigh1;
    xlow2 = onepair->xlow2;
    xhigh2 = onepair->xhigh2;
    xshift = onepair->xshift;

    ylow1 = onepair->ylow1;
    yhigh1 = onepair->yhigh1;
    ylow2 = onepair->ylow2;
    yhigh2 = onepair->yhigh2;
    yshift = onepair->yshift;

    zlow1 = onepair->zlow1;
    zhigh1 = onepair->zhigh1;
    zlow2 = onepair->zlow2;
    zhigh2 = onepair->zhigh2;
    zshift = onepair->zshift;


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

                    *H += inline_ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength1, p3, &ione, p2, &ione);

                    idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                    idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += inline_ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength2, p3, &ione, p2, &ione);

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

                    *H += inline_ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength1, p3, &ione, p2, &ione);

                    idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                    idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += inline_ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength2, p3, &ione, p2, &ione);

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

                    *H += inline_ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength1, p3, &ione, p2, &ione);

                    idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                    idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += inline_ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength2, p3, &ione, p2, &ione);

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

                    *H += inline_ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength1, p3, &ione, p2, &ione);

                    idx1 = ix1 + iy1 + zlow2 - orbit1->izmin - zshift1;
                    idx2 = ix2 + iy2 + zlow2 - orbit2->izmin - zshift2;
                    p1 = &psi1[idx1];
                    p2 = &psi2[idx2];
                    p3 = &psi3[idx1];

                    *H += inline_ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength2, p3, &ione, p2, &ione);

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

                    *H += inline_ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength1, p3, &ione, p2, &ione);

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

                    *H += inline_ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength1, p3, &ione, p2, &ione);

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

                    *H += inline_ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength1, p3, &ione, p2, &ione);

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

                    *H += inline_ddot(&zlength1, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength1, p3, &ione, p2, &ione);

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

                    *H += inline_ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength2, p3, &ione, p2, &ione);

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

                    *H += inline_ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength2, p3, &ione, p2, &ione);

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

                    *H += inline_ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength2, p3, &ione, p2, &ione);

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

                    *H += inline_ddot(&zlength2, p1, &ione, p2, &ione);
                    *S += inline_ddot(&zlength2, p3, &ione, p2, &ione);

                }
            }

            break;

    }

    if(ct.verbose) delete(RT1);

}

