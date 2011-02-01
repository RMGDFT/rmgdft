/************************** SVN Revision Information **************************
 **    $Id: pack_gftoc.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"

void pack_gftoc(SPECIES * sp, fftw_complex * gwptr, fftw_complex * gbptr)
{

    int idx1, idx2, i, j, k, size;
    int i1, i2, j1, j2, k1, k2;

    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    for (i = -sp->nldim / 2; i <= sp->nldim / 2; i++)
    {
        for (j = -sp->nldim / 2; j <= sp->nldim / 2; j++)
        {
            for (k = -sp->nldim / 2; k <= sp->nldim / 2; k++)
            {
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
                gbptr[idx1].re = gwptr[idx2].re / (double) size;
                gbptr[idx1].im = gwptr[idx2].im / (double) size;
            }
        }
    }

}
