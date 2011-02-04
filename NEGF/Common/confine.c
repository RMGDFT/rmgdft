/************************** SVN Revision Information **************************
 **    $Id: confine.c 1242 2011-02-02 18:55:23Z luw $    **
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


void confine (REAL * mat, int size_x, int size_y, int size_z, COMPASS compass, int grid_type)
{
    int i, j, k;
    int x, y, z;
    int pex, pey, pez, xoff, yoff, zoff;
    int test;
    int idx;




    /* find the offset  */
    pe2xyz (pct.thispe, &pex, &pey, &pez);
    xoff = pex * FPX0_GRID;
    yoff = pey * FPY0_GRID;
    zoff = pez * FPZ0_GRID;




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
                    if (grid_type == 2)
                    {
                        x = x - 2;
                        y = y - 2;
                        z = z - 2;
                        if (x < 0)
                            x += FNX_GRID;
                        if (y < 0)
                            y += FNY_GRID;
                        if (z < 0)
                            z += FNZ_GRID;
                        if (x >= FNX_GRID)
                            x -= FNX_GRID;
                        if (y >= FNY_GRID)
                            y -= FNY_GRID;
                        if (z >= FNZ_GRID)
                            z -= FNZ_GRID;
                    }

                    test = ((x < compass.box1.x1) || (x >= compass.box1.x2) ||
                            (y < compass.box1.y1) || (y >= compass.box1.y2) ||
                            (z < compass.box1.z1) || (z >= compass.box1.z2));

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
