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
#include "negf_prototypes.h"
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "twoParts.h"


void confine (double * mat, int size_x, int size_y, int size_z, COMPASS compass, int level)
{
    int i, j, k;
    int x, y, z;
    int pex, pey, pez, xoff, yoff, zoff;
    int test;
    int idx, xgrid, ygrid, zgrid;
    int ioffset, ix, iy,iz;
    double tem;



    /* find the offset  */
    pe2xyz (pct.gridpe, &pex, &pey, &pez);
    xgrid =  get_FNX_GRID()/(1<<level);
    ygrid =  get_FNY_GRID()/(1<<level);
    zgrid =  get_FNZ_GRID()/(1<<level);
    
    xoff = pex * (xgrid /get_PE_X() );
    ix = xgrid % get_PE_X();
    ioffset = 0;

    for(idx = 1;idx <= pex;idx++) {
        if(idx <= ix) ioffset++;
    }
     
    xoff += ioffset;



    yoff = pey * (ygrid /get_PE_Y() );
    iy = ygrid % get_PE_Y();
    ioffset = 0;

    for(idx = 1;idx <= pey;idx++) {
        if(idx <= iy) ioffset++;
    }
     
    yoff += ioffset;

    zoff = pez * (zgrid /get_PE_Z() );
    iz = zgrid % get_PE_Z();
    ioffset = 0;

    for(idx = 1;idx <= pez;idx++) {
        if(idx <= iz) ioffset++;
    }
     
    zoff += ioffset;

    if(level == 0) 
    {
        xoff = get_FPX_OFFSET();
        yoff = get_FPY_OFFSET();
        zoff = get_FPZ_OFFSET();
    }

    if (compass.type == 1)
    {
        for (i = 0; i < size_x; i++)
        {
            x = i + xoff;
            test = ((x <= compass.box1.x1/(1<<level)) || (x >= compass.box1.x2/(1<<level)) );
            if(test)
            {
                for (j = 0; j < size_y; j++)
                {
                    for (k = 0; k < size_z; k++)
                    {
                        idx = i * size_y * size_z + j * size_z + k;

                        mat[idx] = 0.0;

                    }
                }

            }
        }

        for (j = 0; j < size_y; j++)
        {
            y = j + yoff;
            test = ((y <= compass.box1.y1/(1<<level)) || (y >= compass.box1.y2/(1<<level)) );
            if(test)
            {
                for (i = 0; i < size_x; i++)
                {
                    for (k = 0; k < size_z; k++)
                    {
                        idx = i * size_y * size_z + j * size_z + k;

                        mat[idx] = 0.0;

                    }
                }

            }
        }

        for (k = 0; k < size_z; k++)
        {
            z = k + zoff;
            test = ((z <= compass.box1.z1/(1<<level)) || (z >= compass.box1.z2/(1<<level)) );
            if(test)
            {
                for (i = 0; i < size_x; i++)
                {
                    for (j = 0; j < size_y; j++)
                    {
                        idx = i * size_y * size_z + j * size_z + k;

                        mat[idx] = 0.0;

                    }
                }

            }
        }

    }


}
