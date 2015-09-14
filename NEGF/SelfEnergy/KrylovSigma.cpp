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
 *    Krylov method to calculate self energy Sigma
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
#include <complex>
#include "my_scalapack.h"
#include "blas.h"


void KrylovSigma(int n, std::complex<double> *H00, std::complex<double> *H10, 
        std::complex<double> *H01, std::complex<double> *sigma, double lamda_min)
{


    std::complex<double> *A, *W, *VL, *VR, *work;
    std::complex<double> shift, alpha;
    double *rwork;
    int *ipiv, i, j, lwork, info, n2, nn;
    int ione = 1;
    n2 = 2*n;
    nn = n*n;
    lwork = 2* n2;

    A = new std::complex<double>[n2*n2];
    VR = new std::complex<double>[n2*n2];
    W = new std::complex<double>[n2];
    work = new std::complex<double>[lwork];
    rwork = new double[lwork];
    ipiv = new int[n2];

    shift = 1.0/sqrt(2.0);

    // set up matrix A see eq 13 in PRB Kurt

    // first A store M eq 8, VR store C eq.9 and K eq.10
    zcopy(&nn, H10, &ione, A, &ione);
    zaxpy(&nn, &shift, H00, &ione, A, &ione);
    alpha = shift * shift;
    zaxpy(&nn, &alpha, H01, &ione, A, &ione);

    zcopy(&nn, H00, &ione, VR, &ione);

    alpha = 2.0*shift;
    zaxpy(&nn, &alpha, H10, &ione, VR, &ione);
    zcopy(&nn, H10, &ione, &VR[nn], &ione);


    zgesv(&n, &n2, (double *)A, &n, ipiv, (double *)VR, &n, &info );
    if (info != 0)
    {
        printf ("error in zgesv in Krylov with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }

    // now VR store M^-1* C and M^-1 *K  copy them to A
    //  A[j * n2 + i], i row, j col

    for(i = 0; i < n2 * n2; i++) A[i] = 0.0;

    for(i = 0; i < n; i++) A[(n+i) * n2 + i] = 1.0;   

    for(j = 0; j < n; j++) 
        for(i = 0; i < n; i++) 
        {
            A[j * n2 + n+i] = VR[j*n+i];
            A[(j+n) * n2 + n+i] = VR[(j+n) * n + i];
        }


    zgeev("N", "V", &n2, A, &n2, W, VL, &n2, VR, &n2, work, &lwork, rwork, &info );

    if (info != 0)
    {
        printf ("error in zgeev in Krylov with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }

    for (i = 0; i < n2; i++) printf("\n eigv %d %f %f", i, W[i]);
}

