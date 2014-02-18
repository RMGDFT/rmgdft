/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

                            genvpsi.c
 

*/



#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "main_on.h"



void genvnlpsi(double *sg_twovpsipsi, double *vnl, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, ixstart, iystart, ixstartp, iystartp;
    int incx, incy, incx1, incy1;

    incy = dimz + 2;
    incx = (dimy + 2) * (dimz + 2);

    incy1 = dimz;
    incx1 = dimy * dimz;


    /* Generate 2 * V * psi */
    for (ix = 0; ix < dimx; ix++)
    {

        ixstart = ix * incx1;
        ixstartp = (ix + 1) * incx;

        for (iy = 0; iy < dimy; iy++)
        {

            iystart = ixstart + iy * incy1;
            iystartp = ixstartp + (iy + 1) * incy + 1;

            for (iz = 0; iz < dimz; iz++)
            {

                sg_twovpsipsi[iystartp + iz] = 2. * vnl[iystart + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


}                               /* end genvpsi */
