/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <climits>
#include <complex>

void PackGftoc (int nlfxdim, int nlfydim, int nlfzdim, int nlxdim, int nlydim, int nlzdim,
        std::complex<double> *gwptr, std::complex<double> *gbptr)
{

    int idx1, idx2, i, j, k, size;
    int i1, i2, j1, j2, k1, k2;

    int icut = (nlxdim / 2) * (nlxdim / 2);
    
    if(nlxdim != nlydim || nlxdim != nlzdim || nlydim != nlzdim)
    {
        icut = INT_MAX;
        printf("\n WARNING:  not cut in PackGftoc.cpp\n");
    }

    if(nlxdim%2 !=0 || nlydim%2 !=0 || nlzdim%2 !=0)
    {
        printf("\n WARNING:  PackGftoc assum even numbers %d %d %d\n", nlxdim, nlydim, nlzdim);
        fflush(NULL);
        exit(0);
    }
        
        printf("\n WARNING:  PackGftoc assum even numbers %d %d %d\n", nlxdim, nlydim, nlzdim);

    for(idx1 = 0; idx1 < nlxdim * nlydim *nlzdim; idx1++) gbptr[idx1] = 0.0;
    size = nlfxdim * nlfydim * nlfzdim;

    for (i = -nlxdim/2; i < nlxdim/ 2; i++)
    {
        int isq = i * i;
        i1 = (i + nlxdim) %nlxdim;
        i2 = (i + nlfxdim) %nlfxdim;
        for (j = -nlydim/2; j < (nlydim) / 2; j++)
        {
            int jsq = j * j;
            j1 = (j + nlydim) %nlydim;
            j2 = (j + nlfydim) %nlfydim;
            for (k = -nlzdim/2; k < (nlzdim) / 2; k++)
            {
                int ksq = k * k;
                k1 = (k + nlzdim) %nlzdim;
                k2 = (k + nlfzdim) %nlfzdim;


                idx1 = i1 * nlydim * nlzdim + j1 * nlzdim + k1;
                idx2 = i2 * nlfydim * nlfzdim + j2 * nlfzdim + k2;

                if(icut >= (isq + jsq + ksq)) {
                    gbptr[idx1] +=  gwptr[idx2] / (double) size;
                }
            }
        }
    }

}
