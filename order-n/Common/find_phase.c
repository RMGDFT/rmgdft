/************************** SVN Revision Information **************************
 **    $Id: find_phase.c 809 2007-05-29 23:15:00Z miro $    **
******************************************************************************/
 
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "md.h"

/*This calculates phase factor that will be used when calculating backwards fourier transform*/
void find_phase (int nldim, REAL * nlcdrs, REAL * phase_sin, REAL * phase_cos)
{

    int idx1, i, j, k;
    int i1, j1, k1, nldim_sq;
    REAL theta;
    REAL rgs_x, rgs_y, rgs_z;

    /*Reciprocal grid spacings in x, y and z directions */
    rgs_x = 1.0 / (ct.hxgrid * ct.xside);
    rgs_y = 1.0 / (ct.hygrid * ct.yside);
    rgs_z = 1.0 / (ct.hzgrid * ct.zside);

    nldim_sq = nldim * nldim;


    for (i = -nldim / 2; i <= nldim / 2; i++)
    {
        for (j = -nldim / 2; j <= nldim / 2; j++)
        {
            for (k = -nldim / 2; k <= nldim / 2; k++)
            {

                if (i < 0)
                    i1 = i + nldim;
                else
                    i1 = i;

                if (j < 0)
                    j1 = j + nldim;
                else
                    j1 = j;

                if (k < 0)
                    k1 = k + nldim;
                else
                    k1 = k;


                /* Phase factor */
                theta = 2.0 * PI / nldim *
                    (((nlcdrs[0] * (REAL) i) * rgs_x)
                     + ((nlcdrs[1] * (REAL) j) * rgs_y) + ((nlcdrs[2] * (REAL) k) * rgs_z));

                idx1 = i1 * nldim_sq + j1 * nldim + k1;

                phase_sin[idx1] = sin (theta);
                phase_cos[idx1] = cos (theta);
            }
        }
    }

}
