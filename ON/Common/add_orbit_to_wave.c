/************************** SVN Revision Information **************************
 **    $Id: add_orbit_to_wave 2013-06-02 15:38:05Z BTAN $    **
******************************************************************************/
 
/*

Add orbit st1(psi1) multiplied with its coefficient scale to get a particular wave wave_global
wrapper function call in get_wave

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main_on.h"

static int fold_to_unitcell(int, int);

void add_orbit_to_wave(int st1, REAL scale, REAL * psi1, REAL * wave_global, STATE * states)
{

    int iyy, izz, iyy1, izz1;
    int incx, incy;
    int ix, iy, iz;
    int ix1,  iy1,  iz1,  idx1 ;
    int ix3, iy3, iz3, idx3;


    iyy = states[st1].iymax - states[st1].iymin + 1;
    izz = states[st1].izmax - states[st1].izmin + 1;
    incx = iyy * izz;
    incy = izz;


    for (ix = states[st1].ixmin; ix <= states[st1].ixmax; ix++)
    {
        ix1 = (ix - states[st1].ixmin) * incx;//ix1 traverse the relative grids within orbital st1
        ix3 = fold_to_unitcell(ix, get_NX_GRID()) * get_NY_GRID() * get_NZ_GRID();//ix3 traverse the corresponding global coarse grids

        for (iy = states[st1].iymin; iy <= states[st1].iymax; iy++)
        {
            iy1 = (iy - states[st1].iymin) * incy;
            iy3 = fold_to_unitcell(iy, get_NY_GRID()) * get_NZ_GRID();

            for (iz = states[st1].izmin; iz <= states[st1].izmax; iz++)
            {
                iz1 = iz - states[st1].izmin;
                iz3 = fold_to_unitcell(iz, get_NZ_GRID());

                idx1 = ix1 + iy1 + iz1;
                idx3 = ix3 + iy3 + iz3;

                wave_global[idx3] += psi1[idx1]  * scale; // (next time we call this function, we will add another orbital to wave_global)

            }
        }

    }

}



static int fold_to_unitcell(int ix, int NX)
{

    int item;

    if (ix < 0)
        item = ix + NX;
    else if (ix >= NX)
        item = ix - NX;
    else
        item = ix;

    return item;
}
