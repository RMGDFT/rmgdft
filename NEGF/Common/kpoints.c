/*
 **    $Id$    **
******************************************************************************/
 

/*
 *
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <complex.h>



void kpoints(int *nkp, double *kvecx, double *kvecy, double *kvecz, int *nkp_tot1, double *kweight)
{

    int i, j, k, ncheck, kminus;
    double kvec_tem[3], weight, pi;
    int nkp_tot;

    nkp_tot = 0;
    
    pi = 4.0 * atan(1.0);
    for (i = 0; i < nkp[0]; i++)
        for (j = 0; j < nkp[1]; j++)
            for (k = 0; k < nkp[2]; k++)

            {


                kvec_tem[0] =  pi*(1.0 - 1.0/nkp[0] - 2.0 /nkp[0] * i);
                kvec_tem[1] =  pi*(1.0 - 1.0/nkp[1] - 2.0 /nkp[1] * j);
                kvec_tem[2] =  pi*(1.0 - 1.0/nkp[2] - 2.0 /nkp[2] * k);

                // check if this kpoint = -previous kpoint
                
                kminus = 0;
                for (ncheck = 0; ncheck < nkp_tot; ncheck++)
                {
                    if (    fabs(kvec_tem[0] + kvecx[ncheck]) <1.0e-8 &&
                            fabs(kvec_tem[1] + kvecy[ncheck]) <1.0e-8 &&+
                            fabs(kvec_tem[2] + kvecz[ncheck]) <1.0e-8 )
                    {
                        kminus = 1;
                        break;
                    }
                }

                if(kminus) 
                {
                        kweight[ncheck] += 1.0;
                }
                else
                {
                    kvecx[nkp_tot] = kvec_tem[0];
                    kvecy[nkp_tot] = kvec_tem[1];
                    kvecz[nkp_tot] = kvec_tem[2];
                    kweight[nkp_tot] = 1.0;

                    nkp_tot++;

                }

            }

// normalize the weight
    weight = 0.0;
    for(i= 0; i < nkp_tot; i++)
        weight += kweight[i];

    for(i= 0; i < nkp_tot; i++)
        kweight[i] /= weight;

    *nkp_tot1 = nkp_tot;
}

