/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "const.h"
#include "params.h"
#include "rmg_alloc.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"


/*This calculates phase factor that will be used when calculating backwards fourier transform*/
void find_phase (int nldim, double * nlcdrs, double * phase_sin, double * phase_cos)
{

    int idx1, i, j, k;
    int i1, j1, k1, nldim_sq;
    double theta;
    double rgs_x, rgs_y, rgs_z;

    /*Reciprocal grid spacings in x, y and z directions */
    rgs_x = 1.0 / (get_hxgrid() * get_xside());
    rgs_y = 1.0 / (get_hygrid() * get_yside());
    rgs_z = 1.0 / (get_hzgrid() * get_zside());

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
                    (((nlcdrs[0] * (double) i) * rgs_x)
                     + ((nlcdrs[1] * (double) j) * rgs_y) + ((nlcdrs[2] * (double) k) * rgs_z));

                idx1 = i1 * nldim_sq + j1 * nldim + k1;

                phase_sin[idx1] = sin (theta);
                phase_cos[idx1] = cos (theta);
            }
        }
    }

}
