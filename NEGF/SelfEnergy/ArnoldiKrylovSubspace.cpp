/************************** SVN Revision Information **************************
 **    $Id: Main.cpp 3223 2015-09-11 19:05:19Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/md.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 *    generate Krylov subspace with Arnoldi procedure   
 *   
 * INPUTS
 * OUTPUT
 * PARENTS
 * CHILDREN
 * SEE ALSO
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "blas.h"


void ArnoldiKrylovSubspace(int n, int k, int m, std::complex<double> *A, 
        std::complex<double> *V, std::complex<double> *Hkk, std::complex<double> *Hmm)
{

    // expanding the subspace from k vectors to m vectors, m >k
    // A(n,n) matrix, V(n,m), the first V(n,k) is defined  as input
    // Hkk(k,k) matrix, input
    // Hmm(m,m) matrix, its (k,k) block is same as Hkk

    int i,j,seed, ione = 1;
    std::complex<double> norm_coeff, alpha;
    std::complex<double> one(1.0), zero(0.0);
    RmgTimer *RT = new RmgTimer("Generating Subspace");


    // generating the first random vector with length of n
    if(k == 0) 
    {
        seed = 345;
        norm_coeff = 0.0;
        for(i = 0; i <n; i++)
        {
            V[i] = rand0(&seed); 
            norm_coeff += std::norm(V[i]);
        }
        norm_coeff = 1.0/norm_coeff;
        zcal(&n, &norm_coeff; V, &ione);
    }

//   set Hmm = 0 and copy Hkk to Hmm 
    for(i=0; i < m*m; i++) Hmm[i] = 0.0;
    for(i = 0; i < k; i++)
    for(j = 0; j < k; j++)
        Hmm[i*m + j] = Hkk[i*k + j];

    for(j = k; j < m; j++)
    {
        zgemv ("N", &n, &n, &one, A, &n, &V[j*n], &ione, &zero, &V[(j+1)*n], &ione)
            for(i = 0; i <= j; i++)
            {
                Hmm[i*m + j] = zdotc(&n, &V[i*n], &ione, &V[(j+1)*n], &ione);
                alpha = -Hmm[i*m+j];
                zaxpy(&n, &alpha, &V[i*n], &ione, &V[(j+1)*n], &ione);
            }

        if(j < m-1)

        {
            
            Hmm[(j+1) * m + j] = 0.0;
            for(i = 0; i <n; i++)
                Hmm[(j+1) * m + j] += std::norm(V[(j+1)*n + i]);
            norm_coeff = 1.0/ Hmm[(j+1) * m + j];

            zcal(&n, &norm_coeff; &V[(j+1)*m], &ione);
        }
    }
}

