/************************** SVN Revision Information **************************
 **    $Id: Sigma.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/*
 *  Calculate the Sigma matrix,  eq. 20 of PRB 65, 165401   
 *  V = HLC - e * SLC
 *  HLC, SLC:  nL * nC real matrix
 *  Green:  nL * nL complex matrix
 *  sigma: nC * nC complex matrix as output
 *  sigma = V^+ Green * V
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "md.h"


void Sigma (REAL * sigma, REAL * HLC, REAL * SLC, REAL eneR, REAL eneI, REAL * green, int iprobe)
{

    REAL *ch, *ch1, *ch2;
    int i, j, nL, nC, N1;
    char fcd_n = 'N';
    REAL alpha[2], beta[2];


    alpha[0] = 1.0;
    alpha[1] = 0.0;
    beta[0] = 0.0;
    beta[1] = 0.0;

    nL = lcr[iprobe].num_states;
	N1 = cei.probe_in_block[iprobe - 1];
    nC = ct.block_dim[N1];
	
    my_malloc_init( ch, 2 * nL * nC, REAL );
    my_malloc_init( ch1, 2 * nL * nC, REAL );
    my_malloc_init( ch2, 2 * nL * nC, REAL );
    for (i = 0; i < nL * nC; i++)
    {
        ch[i * 2] = HLC[i] * Ha_eV - eneR * SLC[i];
        ch[i * 2 + 1] = -eneI * SLC[i];
    }

    for (i = 0; i < nL; i++)
        for (j = 0; j < nC; j++)
        {
            ch2[i * nC * 2 + 2 * j] = ch[j * nL * 2 + 2 * i];
            ch2[i * nC * 2 + 2 * j + 1] = ch[j * nL * 2 + 2 * i + 1];
        }


/*    ch1(C,L) = ch^+(C,L) * green(L,L)   */
/*    sigma(C,C) = ch1(C,L) * ch(L,C)   */
    
        QMD_ZGEMM (&fcd_n, &fcd_n, &nC, &nL, &nL, alpha, ch, &nC, green, &nL, beta, ch1, &nC);

        QMD_ZGEMM (&fcd_n, &fcd_n, &nC, &nC, &nL, alpha, ch1, &nC, ch2, &nL, beta, sigma, &nC);
   

    my_free(ch);
    my_free(ch1);
    my_free(ch2);

}
