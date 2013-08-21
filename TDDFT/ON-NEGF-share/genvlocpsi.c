/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

                            genvlocpsi.c
 

*/



#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "main.h"



void genvlocpsi(REAL * psi, int st1, REAL * work1, REAL * vtot_global, STATE * states)
{

    int ix, iy, iz, ixstart, iystart, ixstartp, iystartp;
    int incx, incy;
    int incx1, incy1;
    int ix1, iy1, iz1, ix2, iy2, iz2;
    int idx1, idx2;

    ix1 = states[st1].ixmin;
    iy1 = states[st1].iymin;
    iz1 = states[st1].izmin;
    ix2 = states[st1].ixmax;
    iy2 = states[st1].iymax;
    iz2 = states[st1].izmax;

    incx = (iy2 - iy1 + 1) * (iz2 - iz1 + 1);
    incy = (iz2 - iz1 + 1);

    incx1 = NY_GRID * NZ_GRID;
    incy1 = NZ_GRID;


    /* Generate 2 * V * psi */
    for (ix = ix1; ix <= ix2; ix++)
    {

        ixstart = (ix - ix1) * incx;
        if (ix < 0)
        {
            ixstartp = (ix + NX_GRID) * incx1;
        }
        else if (ix >= NX_GRID)
        {
            ixstartp = (ix - NX_GRID) * incx1;
        }
        else
        {
            ixstartp = ix * incx1;
        }

        for (iy = iy1; iy <= iy2; iy++)
        {

            iystart = ixstart + (iy - iy1) * incy;
            if (iy < 0)
            {
                iystartp = ixstartp + (iy + NY_GRID) * incy1;
            }
            else if (iy >= NY_GRID)
            {
                iystartp = ixstartp + (iy - NY_GRID) * incy1;
            }
            else
            {
                iystartp = ixstartp + iy * incy1;
            }

            for (iz = iz1; iz <= iz2; iz++)
            {
                idx1 = iystart + iz - iz1;
                if (iz < 0)
                {
                    idx2 = iystartp + (iz + NZ_GRID);
                }
                else if (iz >= NZ_GRID)
                {
                    idx2 = iystartp + (iz - NZ_GRID);
                }
                else
                {
                    idx2 = iystartp + iz;
                }

                work1[idx1] = 2.0 * psi[idx1] * vtot_global[idx2];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


}                               /* end genvlocpsi */
