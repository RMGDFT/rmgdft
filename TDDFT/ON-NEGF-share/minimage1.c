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

REAL minimage1(REAL aa[3], REAL bb[3])
{

    int ix, iy, iz, idx, idxmin, nn = 0;
    REAL r[27], ax, ay, az, rmin, x, y, z;


    /* Get lattice vectors */
    ax = NX_GRID * ct.hxgrid * ct.xside;
    ay = NY_GRID * ct.hygrid * ct.yside;
    az = NZ_GRID * ct.hzgrid * ct.zside;



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

                    x = aa[0] - bb[0] + (REAL) ix *ax;
                    y = aa[1] - bb[1] + (REAL) iy *ay;
                    z = aa[2] - bb[2] + (REAL) iz *az;

                    r[idx] = (x * x + y * y + z * z);

                    idx++;

                }               /* end for */

            }                   /* end for */

        }                       /* end for */
    }
    else if (ct.boundaryflag == ZPERIODIC)
    {
        nn = 3;

        idx = 0;
        for (iz = -1; iz <= 1; iz++)
        {

            x = aa[0] - bb[0];
            y = aa[1] - bb[1];
            z = aa[2] - bb[2] + (REAL) iz *az;

            r[idx] = (x * x + y * y + z * z);

            idx++;

        }                       /* end for */
    }


    /* Next we find rmin */
    rmin = 10000000.0;
    for (idx = 0; idx < nn; idx++)
    {

        if (r[idx] < rmin)
        {

            rmin = r[idx];
            idxmin = idx;

        }                       /* end if */

    }                           /* end for */

    rmin = sqrt(rmin);
    return rmin;

}                               /* end minimage1 */
