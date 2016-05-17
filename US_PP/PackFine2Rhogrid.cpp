/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

void PackFine2Rhogrid (std::complex<double> *gwptr, int ngrid_fine, std::complex<double> *gbptr, int ngrid) 
{

    int idx1, idx2, i, j, k, size;
    int i1, i2, j1, j2, k1, k2;

    size = ngrid_fine * ngrid_fine * ngrid_fine;
    int icut = (ngrid / 2) * (ngrid/ 2);

    for (i = -ngrid / 2; i <= ngrid / 2; i++)
    {
        int isq = i * i;
        for (j = -ngrid / 2; j <= ngrid / 2; j++)
        {
            int jsq = j * j;
            for (k = -ngrid / 2; k <= ngrid / 2; k++)
            {
                int ksq = k * k;
                if (i < 0)
                {
                    i1 = i + ngrid;
                    i2 = i + ngrid_fine;
                }
                else
                {
                    i1 = i;
                    i2 = i;
                }
                if (j < 0)
                {
                    j1 = j + ngrid;
                    j2 = j + ngrid_fine;
                }
                else
                {
                    j1 = j;
                    j2 = j;
                }
                if (k < 0)
                {
                    k1 = k + ngrid;
                    k2 = k + ngrid_fine;
                }
                else
                {
                    k1 = k;
                    k2 = k;
                }
                idx1 = i1 * ngrid * ngrid + j1 * ngrid + k1;
                idx2 = i2 * ngrid_fine * ngrid_fine + j2 * ngrid_fine + k2;
                
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
