/************************** SVN Revision Information **************************
 **    $Id: diff_hx_interpolation.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/*
 *  interpolation routine for rho, potentials vh and vxc
 *  in x direction, number of grids is same, but the hx may be slightly 
 *  different after we conbine the lead and conductor togather 
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "md.h"

void diff_hx_interpolation (double *xi, double *xi_old, int NX,
                            double hx, double hx_old, double x0, double x0_old)
{

    int ix, i1, i0, i2, i3;
    double x, x1, frac;
    double cc0, cc1, cc2, cc3;


    for (ix = 0; ix < NX; ix++)
    {
        x = x0 + ix * hx;

        x1 = (x - x0_old) / hx_old;
        i1 = (int) x1;

        if (x1 < -2.0)
        {
            printf ("\n x1 = %f %f %f %f %f", x1, x0, x0_old, hx, hx_old);
            error_handler ("grid ???");
        }


        if (x1 < 0.0) i1 = -1;
        if (x1 < -1.0) i1 = -2;

        frac = (x - i1 * hx_old - x0_old) / hx_old;

        i0 = i1 - 1;
        i2 = i1 + 1;
        i3 = i1 + 2;

/* fold if the points are around boarder */

        if (i0 < 0)
            i0 += NX;
        if (i1 < 0)
            i1 += NX;
        if (i2 < 0)
            i2 += NX;
        if (i3 < 0)
            i3 += NX;
        if (i0 >= NX)
            i0 -= NX;
        if (i1 >= NX)
            i1 -= NX;
        if (i2 >= NX)
            i2 -= NX;
        if (i3 >= NX)
            i3 -= NX;

        cc0 = -frac * (1.0 - frac) * (2.0 - frac) / 6.0;
        cc1 = (1.0 + frac) * (1.0 - frac) * (2.0 - frac) / 2.0;
        cc2 = (1.0 + frac) * frac * (2.0 - frac) / 2.0;
        cc3 = -(1.0 + frac) * frac * (1.0 - frac) / 6.0;

        xi[ix] = cc0 * xi_old[i0] + cc1 * xi_old[i1] + cc2 * xi_old[i2] + cc3 * xi_old[i3];
/*	if(pct.thispe ==0) printf("\n %d %f %f dddd", ix, xi[ix], xi_old[ix]); */
    }


}
