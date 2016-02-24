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

    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;
    int icut = (sp->nldim / 2) * (sp->nldim / 2);

    for (i = -sp->nldim / 2; i <= sp->nldim / 2; i++)
    {
        int isq = i * i;
        for (j = -sp->nldim / 2; j <= sp->nldim / 2; j++)
        {
            int jsq = j * j;
            for (k = -sp->nldim / 2; k <= sp->nldim / 2; k++)
            {
                int ksq = k * k;
                if (i < 0)
                {
                    i1 = i + sp->nldim;
                    i2 = i + sp->nlfdim;
                }
                else
                {
                    i1 = i;
                    i2 = i;
                }
                if (j < 0)
                {
                    j1 = j + sp->nldim;
                    j2 = j + sp->nlfdim;
                }
                else
                {
                    j1 = j;
                    j2 = j;
                }
                if (k < 0)
                {
                    k1 = k + sp->nldim;
                    k2 = k + sp->nlfdim;
                }
                else
                {
                    k1 = k;
                    k2 = k;
                }
                idx1 = i1 * sp->nldim * sp->nldim + j1 * sp->nldim + k1;
                idx2 = i2 * sp->nlfdim * sp->nlfdim + j2 * sp->nlfdim + k2;
                
                if(icut > (isq + jsq + ksq)) {
                    gbptr[idx1] =  gwptr[idx2] / (double) size;
                }
                else {
                    gbptr[idx1] = 0.0 + 0.0I;
                }
            }
        }
    }

}
