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

    int nlxdim = sp->nldim;
    int nlydim = sp->nldim;
    int nlzdim = sp->nldim;
    if(!ct.localize_projectors) {
        nlxdim = get_NX_GRID();
        nlydim = get_NY_GRID();
        nlzdim = get_NZ_GRID();
    }

    idx = nlxdim * nlydim * nlzdim;
    my_malloc(tem_weptr, idx, fftw_complex);
    
    for(ix = 0; ix < idx; ix++) tem_weptr[ix] = weptr[ix];

    for (ix = 0; ix < nlxdim; ix++)
    {
        ixx = ix + nlxdim/2;
        if(ixx >= nlxdim) ixx = ixx - nlxdim;

        for (iy = 0; iy < nlydim; iy++)
        {
            iyy = iy + nlydim/2;
            if(iyy >= nlydim) iyy = iyy - nlydim;

            for (iz = 0; iz < nlzdim; iz++)
            {
                izz = iz + nlzdim/2;
                if(izz >= nlzdim) izz = izz - nlzdim;

                idx = ix * nlydim * nlzdim + iy * nlzdim + iz;
                idx1 = ixx * nlydim * nlzdim + iyy * nlzdim + izz;
                weptr[idx1] = tem_weptr[idx];

            }                   /* end for */
        }                       /* end for */
    }                           /* end for */

    my_free(tem_weptr);

}
