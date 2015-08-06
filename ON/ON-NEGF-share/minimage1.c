/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
                            make_mask_grid.c
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"

double minimage1(double aa[3], double bb[3])
{

    int ix, iy, iz, idx, nn = 0;
    double r[27], ax, ay, az, rmin, x, y, z;


    /* Get lattice vectors */
    ax = get_NX_GRID() * get_hxgrid() * get_xside();
    ay = get_NY_GRID() * get_hygrid() * get_yside();
    az = get_NZ_GRID() * get_hzgrid() * get_zside();



    /* Loop over all possible combinations and generate r */
    if (ct.boundaryflag == PERIODIC)
    {
        nn = 27;

        idx = 0;
        for (ix = -1; ix <= 1; ix++)
        {

            for (iy = -1; iy <= 1; iy++)
            {

                for (iz = -1; iz <= 1; iz++)
                {

                    x = aa[0] - bb[0] + (double) ix *ax;
                    y = aa[1] - bb[1] + (double) iy *ay;
                    z = aa[2] - bb[2] + (double) iz *az;

                    r[idx] = (x * x + y * y + z * z);

                    idx++;

                }               /* end for */

            }                   /* end for */

        }                       /* end for */
    }


    /* Next we find rmin */
    rmin = 10000000.0;
    for (idx = 0; idx < nn; idx++)
    {

        if (r[idx] < rmin)
        {

            rmin = r[idx];

        }                       /* end if */

    }                           /* end for */

    rmin = sqrt(rmin);
    return rmin;

}                               /* end minimage1 */
