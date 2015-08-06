/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "main.h"
#include "prototypes_on.h"

void weight_shift(SPECIES *sp, fftw_complex *weptr);
void assign_weight_on(SPECIES * sp, fftw_complex * weptr, double * rtptr)
{

    int idx, ix, iy, iz;

    weight_shift_center(sp, weptr);
    idx = 0;
    for (ix = 0; ix < sp->nldim; ix++)
    {

        for (iy = 0; iy < sp->nldim; iy++)
        {

            for (iz = 0; iz < sp->nldim; iz++)
            {

                rtptr[idx] = creal(weptr[idx]);
//                printf("\n dddddd %d %d %d %f", ix, iy, iz, rtptr[idx]); 
                if (fabs(cimag(weptr[idx])) > 1.0e-6)
                {
                    printf("weptr[%d].im=%e\n", idx, cimag(weptr[idx]));
                    error_handler("something wrong with the fourier transformation");
                }

                idx++;
            }                   /* end for */
        }                       /* end for */
    }                           /* end for */
}
