/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"

/*This is almost exact copy of assign_weight, except that input is real array, not fftw_complex*/

void assign_weight2 (int nldim, int ion, rmg_double_t * beptr, rmg_double_t * rtptr)
{
#if FDIFF_BETA

    int idx, ix, iy, iz, *dvec;
    int idx1, docount;
    int inc_x, inc_y;

    inc_x = nldim * nldim;
    inc_y = nldim;

    dvec = pct.idxflag[ion];
    idx = docount = 0;
    for (ix = 0; ix < nldim; ix++)
    {

        for (iy = 0; iy < nldim; iy++)
        {

            for (iz = 0; iz < nldim; iz++)
            {

                if (dvec[idx])
                {
                    idx1 = ix * inc_x + iy * inc_y + iz;
                    rtptr[docount] = beptr[idx1];
                    docount++;
                }

                idx++;
            }                   /* end for */
        }                       /* end for */
    }                           /* end for */
    if (docount != pct.idxptrlen[ion])
        error_handler ("wrong numbers of projectors");
#endif
}
