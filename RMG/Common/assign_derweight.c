/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"

void assign_derweight (SPECIES * sp, int ion, fftw_complex * beptr, rmg_double_t * rtptr)
{

    int idx, ix, iy, iz, *dvec;
    int idx1, docount;

    dvec = pct.idxflag[ion];
    idx = docount = 0;
    for (ix = 0; ix < sp->nldim; ix++)
    {

        for (iy = 0; iy < sp->nldim; iy++)
        {

            for (iz = 0; iz < sp->nldim; iz++)
            {

                if (dvec[idx])
                {
                    idx1 = ix * sp->nldim * sp->nldim + iy * sp->nldim + iz;
                    rtptr[docount] = beptr[idx1].re;
                    if (beptr[idx1].im > 1.0e-8)
                    {
                        printf ("beptr[%d].im=%e\n", idx1, beptr[idx1].im);
                        error_handler ("something wrong with the fourier transformation");
                    }
                    docount++;
                }

                idx++;
            }
        }
    }
    if (docount != pct.idxptrlen[ion])
    {
        printf ("docount = %d != %d = pct.idxptrlen[ion = %d]\n", docount, pct.idxptrlen[ion], ion);
        error_handler ("wrong numbers of projectors");
    }
}
