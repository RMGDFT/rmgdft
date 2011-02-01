/************************** SVN Revision Information **************************
 **    $Id: assign_weight.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"

void assign_weight(SPECIES * sp, fftw_complex * weptr, REAL * rtptr)
{

    int idx, ix, iy, iz;

    idx = 0;
    for (ix = 0; ix < sp->nldim; ix++)
    {

        for (iy = 0; iy < sp->nldim; iy++)
        {

            for (iz = 0; iz < sp->nldim; iz++)
            {

                rtptr[idx] = weptr[idx].re;
                if (fabs(weptr[idx].im) > 1.0e-6)
                {
                    printf("weptr[%d].im=%e\n", idx, weptr[idx].im);
                    error_handler("something wrong with the fourier transformation");
                }

                idx++;
            }                   /* end for */
        }                       /* end for */
    }                           /* end for */
}
