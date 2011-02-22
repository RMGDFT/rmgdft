/************************** SVN Revision Information **************************
 **    $Id: state_corner_xyz.c 1059 2009-07-27 19:10:28Z luw $    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void state_corner_xyz(STATE * states)
{
    int ix, iy, iz, state;
    REAL hgrid[3];
    REAL tem;
    int ixx, iyy, izz, imin, imax;
    int item;
    int max_nx, max_ny, max_nz;

    hgrid[0] = ct.hxgrid * ct.xside;
    hgrid[1] = ct.hygrid * ct.yside;
    hgrid[2] = ct.hzgrid * ct.zside;



    /* Loop over states */
    for (state = 0; state < ct.num_states; state++)
    {

        /* find the minmum and maximum grid point of the orbit on ion */
        /* x direction */

        imin = 100000;
        imax = -100000;
        tem = states[state].crds[0];

        for (ix = -NX_GRID; ix < 2 * NX_GRID - 1; ix++)
        {
            if (fabs(tem - hgrid[0] * ix) < states[state].radius)
            {
                imin = min(imin, ix);
                imax = max(imax, ix);
            }
        }

        states[state].ixmin = imin - 1;
        item = imax - imin + 1;
        if (item / 2 * 2 < item)
            imax += 1;
        states[state].ixmax = imax + 1;

        /*  y direction */

        imin = 100000;
        imax = -100000;
        tem = states[state].crds[1];

        for (iy = -NY_GRID; iy < 2 * NY_GRID - 1; iy++)
        {
            if (fabs(tem - hgrid[1] * iy) < states[state].radius)
            {
                imin = min(imin, iy);
                imax = max(imax, iy);
            }                   /* endif */
        }                       /* end for */

        states[state].iymin = imin - 1;
        item = imax - imin + 1;
        if (item / 2 * 2 < item)
            imax += 1;
        states[state].iymax = imax + 1;

        /*  z direction */

        imin = 100000;
        imax = -100000;
        tem = states[state].crds[2];

        for (iz = -NZ_GRID; iz < 2 * NZ_GRID - 1; iz++)
        {
            if (fabs(tem - hgrid[2] * iz) < states[state].radius)
            {
                imin = min(imin, iz);
                imax = max(imax, iz);
            }
        }

        states[state].izmin = imin - 1;
        item = imax - imin + 1;
        if (item / 2 * 2 < item)
            imax += 1;
        states[state].izmax = imax + 1;
    }

    for (state = 0; state < ct.num_states; state++)
    {
        ixx = states[state].ixmax - states[state].ixmin + 1;
        iyy = states[state].iymax - states[state].iymin + 1;
        izz = states[state].izmax - states[state].izmin + 1;
        if (ixx > states[state].orbit_nx || iyy > states[state].orbit_ny
                || izz > states[state].orbit_nz)
        {
            printf
                ("states[state].orbit_nx, ny, nz and ixx: %d, %d, %d,   %d \n",
                 states[state].orbit_nx, states[state].orbit_ny, states[state].orbit_nz, ixx);
            error_handler("ixx > states[state].orbit_nx/ny/nz");
        }

        if (states[state].orbit_nx - ixx == 1)
            states[state].ixmin -= 1;
        if (states[state].orbit_nx - ixx == 2)
        {
            states[state].ixmin -= 1;
            states[state].ixmax += 1;
        }
        if (states[state].orbit_nx - ixx == 3)
        {
            states[state].ixmin -= 2;
            states[state].ixmax += 1;
        }
        if (states[state].orbit_nx - ixx == 4)
        {
            states[state].ixmin -= 2;
            states[state].ixmax += 2;
        }
        if (states[state].orbit_nx - ixx > 4)
            error_handler(" states[state].orbit_nx - ixx >4");

        if (states[state].orbit_ny - iyy == 1)
            states[state].iymin -= 1;
        if (states[state].orbit_ny - iyy == 2)
        {
            states[state].iymin -= 1;
            states[state].iymax += 1;
        }
        if (states[state].orbit_ny - iyy == 3)
        {
            states[state].iymin -= 2;
            states[state].iymax += 1;
        }
        if (states[state].orbit_ny - iyy == 4)
        {
            states[state].iymin -= 2;
            states[state].iymax += 2;
        }
        if (states[state].orbit_ny - iyy > 4)
            error_handler(" states[state].orbit_ny - iyy >4");

        if (states[state].orbit_nz - izz == 1)
            states[state].izmin -= 1;
        if (states[state].orbit_nz - izz == 2)
        {
            states[state].izmin -= 1;
            states[state].izmax += 1;
        }
        if (states[state].orbit_nz - izz == 3)
        {
            states[state].izmin -= 2;
            states[state].izmax += 1;
        }
        if (states[state].orbit_nz - izz == 4)
        {
            states[state].izmin -= 2;
            states[state].izmax += 2;
        }
        if (states[state].orbit_nz - izz > 4)
            error_handler(" states[state].orbit_nz - izz >4");
    }


    /* include orbit size in x,y,z by 8 */

    int maxnx = 0; 
    int maxny = 0;
    int maxnz = 0;
    int minnx = 1000000; 
    int minny = 1000000;
    int minnz = 1000000;
    for (state = 0; state < ct.num_states; state++)
    {

        states[state].ixmin -= 2; 
        states[state].ixmax += 2; 
        states[state].iymin -= 2; 
        states[state].iymax += 2; 
        states[state].izmin -= 2; 
        states[state].izmax += 2; 

        states[state].orbit_nx +=4; 
        states[state].orbit_ny +=4; 
        states[state].orbit_nz +=4; 

        states[state].size = states[state].orbit_nx * states[state].orbit_ny * states[state].orbit_ny;    
        maxnx = max(maxnx, states[state].orbit_nx);
        maxny = max(maxny, states[state].orbit_ny);
        maxnz = max(maxnz, states[state].orbit_nz);
        minnx = min(minnx, states[state].orbit_nx);
        minny = min(minny, states[state].orbit_ny);
        minnz = min(minnz, states[state].orbit_nz);

    }

    if(pct.thispe == 0)
    {
        printf("\n lagest orbital size: %d %d %d", maxnx, maxny, maxnz);
        printf("\n smallest orbital size: %d %d %d", minnx, minny, minnz);
    }

    pct.psi_size = 0;
    max_nx = 0;
    max_ny = 0;
    max_nz = 0;
    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        max_nx = max(max_nx, states[state].orbit_nx);
        max_ny = max(max_ny, states[state].orbit_ny);
        max_nz = max(max_nz, states[state].orbit_nz);
        pct.psi_size += states[state].size;
    }
    ct.max_orbit_size = max_nx * max_ny * max_nz +1;
    ct.max_orbit_nx = max_nx;
    ct.max_orbit_ny = max_ny;
    ct.max_orbit_nz = max_nz;

    for (state = 0; state < ct.num_states; state++)
    {
        if (ct.xside <= 2.0 * states[state].radius)
        {
            states[state].ixmin = 0;
            states[state].ixmax = NX_GRID - 1;
            states[state].xfold = 1;
        }
    }
    for (state = 0; state < ct.num_states; state++)
    {
        if (ct.yside <= 2.0 * states[state].radius)
        {
            states[state].iymin = 0;
            states[state].iymax = NY_GRID - 1;
            states[state].yfold = 1;
        }
    }
    for (state = 0; state < ct.num_states; state++)
    {
        if (ct.zside <= 2.0 * states[state].radius)
        {
            states[state].izmin = 0;
            states[state].izmax = NZ_GRID - 1;
            states[state].zfold = 1;
        }
    }
}
