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
    if(!ct.localize_projectors) {
        nlfxdim = ct.nxfgrid * get_NX_GRID();
        nlfydim = ct.nxfgrid * get_NY_GRID();
        nlfzdim = ct.nxfgrid * get_NZ_GRID();
        nlxdim = get_NX_GRID();
        nlydim = get_NY_GRID();
        nlzdim = get_NZ_GRID();
    }

    size = nlfxdim * nlfydim * nlfzdim;
    int icut = (sp->nldim / 2) * (sp->nldim / 2);

    for (i = -(nlxdim-1) / 2; i <= (nlxdim-1) / 2; i++)
    {
        int isq = i * i;
        for (j = -(nlydim-1) / 2; j <= (nlydim-1) / 2; j++)
        {
            int jsq = j * j;
            for (k = -(nlzdim-1) / 2; k <= (nlzdim-1) / 2; k++)
            {
                int ksq = k * k;
                if (i < 0)
                {
                    i1 = i + nlxdim;
                    i2 = i + nlfxdim;
                }
                else
                {
                    i1 = i;
                    i2 = i;
                }
                if (j < 0)
                {
                    j1 = j + nlydim;
                    j2 = j + nlfydim;
                }
                else
                {
                    j1 = j;
                    j2 = j;
                }
                if (k < 0)
                {
                    k1 = k + nlzdim;
                    k2 = k + nlfzdim;
                }
                else
                {
                    k1 = k;
                    k2 = k;
                }
                idx1 = i1 * nlydim * nlzdim + j1 * nlzdim + k1;
                idx2 = i2 * nlfydim * nlfzdim + j2 * nlfzdim + k2;
                
                if(icut >= (isq + jsq + ksq)) {
                    gbptr[idx1] =  gwptr[idx2] / (double) size;
                }
                else {
                    gbptr[idx1] = 0.0 + 0.0I;
                }
            }
        }
    }

}
