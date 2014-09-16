/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "main.h"

void weight_shift_center(SPECIES * sp, fftw_complex * weptr)
{

    int idx, idx1,  ix, iy, iz;
    int ixx, iyy, izz;
    fftw_complex *tem_weptr;

    idx = sp->nldim * sp->nldim * sp->nldim;
    my_malloc(tem_weptr, idx, fftw_complex);
    
    for(ix = 0; ix < idx; ix++) tem_weptr[ix] = weptr[ix];

    for (ix = 0; ix < sp->nldim; ix++)
    {
        ixx = ix + sp->nldim/2;
        if(ixx >=sp->nldim) ixx = ixx - sp->nldim;

        for (iy = 0; iy < sp->nldim; iy++)
        {
            iyy = iy + sp->nldim/2;
            if(iyy >=sp->nldim) iyy = iyy-  sp->nldim;

            for (iz = 0; iz < sp->nldim; iz++)
            {
                izz = iz + sp->nldim/2;
                if(izz >=sp->nldim) izz = izz- sp->nldim;

                idx = ix * sp->nldim * sp->nldim + iy * sp->nldim + iz;
                idx1 = ixx * sp->nldim * sp->nldim + iyy * sp->nldim + izz;
                weptr[idx1] = tem_weptr[idx];

            }                   /* end for */
        }                       /* end for */
    }                           /* end for */

    my_free(tem_weptr);

}
