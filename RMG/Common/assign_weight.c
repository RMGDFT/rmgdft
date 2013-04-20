/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"

void assign_weight (SPECIES * sp, int ion, fftw_complex * beptr, REAL * rtptr)
{

    int idx, ix, iy, iz, *dvec;
    int idx1, docount;
    int *pidx;

    for(idx = 0; idx < pct.P0_BASIS; idx++) rtptr[idx] = 0.0;
    if(pct.idxptrlen[ion] == 0) return;
    pidx = pct.nlindex[ion];
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
                    rtptr[pidx[docount]] = beptr[idx1].re;
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
