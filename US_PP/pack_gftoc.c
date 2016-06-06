/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "main.h"

void pack_gftoc (SPECIES * sp, double complex * gwptr, double complex * gbptr)
{

    int idx1, idx2, i, j, k, size;
    int i1, i2, j1, j2, k1, k2;

    int nlfxdim = sp->nlfdim;
    int nlfydim = sp->nlfdim;
    int nlfzdim = sp->nlfdim;
    int nlxdim = sp->nldim;
    int nlydim = sp->nldim;
    int nlzdim = sp->nldim;
    int icut = (sp->nldim / 2) * (sp->nldim / 2);
    int ixstart = -(nlxdim )/2;
    int iystart = -(nlydim )/2;
    int izstart = -(nlzdim )/2;

    if(!ct.localize_projectors) {
        nlfxdim = ct.nxfgrid * get_NX_GRID();
        nlfydim = ct.nxfgrid * get_NY_GRID();
        nlfzdim = ct.nxfgrid * get_NZ_GRID();
        nlxdim = get_NX_GRID();
        nlydim = get_NY_GRID();
        nlzdim = get_NZ_GRID();
        int mindim = nlxdim;
        if(nlydim < mindim) mindim = nlydim;
        if(nlzdim < mindim) mindim = nlzdim;
        icut = (mindim / 2 - 1) * (mindim / 2 - 1);

        ixstart = -nlxdim/2;
        iystart = -nlydim/2;
        izstart = -nlzdim/2;
    
    }

    for(idx1 = 0; idx1 < nlxdim * nlydim *nlzdim; idx1++) gbptr[idx1] = 0.0 + 0.0I;
    size = nlfxdim * nlfydim * nlfzdim;

    for (i = -nlxdim/2; i < nlxdim/ 2; i++)
    {
        int isq = i * i;
        i1 = (i + nlxdim) %nlxdim;
        i2 = (i + nlfxdim) %nlfxdim;
        for (j = iystart; j < (nlydim) / 2; j++)
        {
            int jsq = j * j;
            j1 = (j + nlydim) %nlydim;
            j2 = (j + nlfydim) %nlfydim;
            for (k = izstart; k < (nlzdim) / 2; k++)
            {
                int ksq = k * k;
                k1 = (k + nlzdim) %nlzdim;
                k2 = (k + nlfzdim) %nlfzdim;


                idx1 = i1 * nlydim * nlzdim + j1 * nlzdim + k1;
                idx2 = i2 * nlfydim * nlfzdim + j2 * nlfzdim + k2;

                //if(icut >= (isq + jsq + ksq)) {
                {
                    //gbptr[idx1] =  wx*wy*wz*gwptr[idx2] / (double) size;
                    gbptr[idx1] +=  gwptr[idx2] / (double) size;
                }
            }
        }
    }

}
