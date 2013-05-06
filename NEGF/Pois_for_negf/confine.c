/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
****f* QMD-MGDFT/confine.c *****
 * 
 */



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"


void confine (REAL * mat, int size_x, int size_y, int size_z, COMPASS compass, int level)
{
    int i, j, k;
    int x, y, z;
    int pex, pey, pez, xoff, yoff, zoff;
    int test;
    int idx, xgrid, ygrid, zgrid;
    int ioffset, ix, iy,iz;
    double tem;



    double time1, time2;

    time1 = my_crtc();
    /* find the offset  */
    pe2xyz (pct.gridpe, &pex, &pey, &pez);
    xgrid =  FNX_GRID/(1<<level);
    ygrid =  FNY_GRID/(1<<level);
    zgrid =  FNZ_GRID/(1<<level);
    
    xoff = pex * (xgrid /PE_X );
    ix = xgrid % PE_X;
    ioffset = 0;

    for(idx = 1;idx <= pex;idx++) {
        if(idx <= ix) ioffset++;
    }
     
    xoff += ioffset;

    yoff = pey * (ygrid /PE_Y );
    iy = ygrid % PE_Y;
    ioffset = 0;

    for(idx = 1;idx <= pey;idx++) {
        if(idx <= iy) ioffset++;
    }
     
    yoff += ioffset;

    zoff = pez * (zgrid /PE_Z );
    iz = zgrid % PE_Z;
    ioffset = 0;

    for(idx = 1;idx <= pez;idx++) {
        if(idx <= iz) ioffset++;
    }
     
    zoff += ioffset;



    if (compass.type == 1)
        for (i = 0; i < size_x; i++)
        {
            for (j = 0; j < size_y; j++)
            {
                for (k = 0; k < size_z; k++)
                {
                    idx = i * size_y * size_z + j * size_z + k;
                    x = i + xoff;
                    y = j + yoff;
                    z = k + zoff;

                    test = ((x <= compass.box1.x1/(1<<level)) || (x >= compass.box1.x2/(1<<level)) ||
                            (y <= compass.box1.y1/(1<<level)) || (y >= compass.box1.y2/(1<<level)) ||
                            (z <= compass.box1.z1/(1<<level)) || (z >= compass.box1.z2/(1<<level)));

                    if (!test)
                    {
                    }
                    else
                    {
                        mat[idx] = 0.0;
                    }


                }
            }

        }

    time2 = my_crtc();

    rmg_timings (CONFINE_TIME, time2 - time1);

}
