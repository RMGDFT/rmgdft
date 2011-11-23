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
#include "md.h"


void confine (REAL * mat, int size_x, int size_y, int size_z, COMPASS compass, int level)
{
    int i, j, k;
    int x, y, z;
    int pex, pey, pez, xoff, yoff, zoff;
    int test;
    int idx;
    double tem;




    /* find the offset  */
    pe2xyz (pct.gridpe, &pex, &pey, &pez);
    xoff = pex * FPX0_GRID/(1<<level);
    yoff = pey * FPY0_GRID/(1<<level);
    zoff = pez * FPZ0_GRID/(1<<level);




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

                    test = ((x < compass.box1.x1/(1<<level)) || (x >= compass.box1.x2/(1<<level)) ||
                            (y < compass.box1.y1/(1<<level)) || (y >= compass.box1.y2/(1<<level)) ||
                            (z < compass.box1.z1/(1<<level)) || (z >= compass.box1.z2/(1<<level)));

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

}
