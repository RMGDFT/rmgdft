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


#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "LCR.h"
#include "prototypes_on.h"
#include "prototypes_negf.h"
#include "init_var.h"



#include "Scalapack.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"


extern "C" void KrylovSigma_c(int n, std::complex<double> *H00, std::complex<double> *H10, 
        std::complex<double> *H01, std::complex<double> *sigma, double lamda_min)
{
    KrylovSigma(n, H00, H10, H01, sigma, lamda_min);
}

void KrylovSigma(int n, std::complex<double> *H00, std::complex<double> *H10, 
        std::complex<double> *H01, std::complex<double> *sigma, double lamda_min)
{


    std::complex<double> *A, *W, *VL, *VR, *work, *lamda;
    std::complex<double> shift, alpha, zero=0.0;
    double *rwork, tol = 1.0e-5;
    int *ipiv, i, j, lwork, info, n2, nn;
    int ione = 1, num_modes;
    n2 = 2*n;
    nn = n*n;
    lwork = 2* n2;

    A = new std::complex<double>[n2*n2];
    VR = new std::complex<double>[n2*n2];
    W = new std::complex<double>[n2];
    lamda = new std::complex<double>[n2];
    work = new std::complex<double>[lwork];
    rwork = new double[lwork];
    ipiv = new int[n2];

    shift.imag( 1.0/sqrt(2.0));
    // set up matrix A see eq 13 in PRB Kurt

    // first A store M eq 8, VR store C eq.9 and K eq.10
    zcopy(&nn, H10, &ione, A, &ione);
    zaxpy(&nn, &shift, H00, &ione, A, &ione);
    alpha = shift * shift;
    zaxpy(&nn, &alpha, H01, &ione, A, &ione);

    zcopy(&nn, H00, &ione, &VR[nn], &ione);

    alpha = 2.0*shift;
    zaxpy(&nn, &alpha, H01, &ione, &VR[nn], &ione);
    zcopy(&nn, H01, &ione, VR, &ione);


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
            A[j * n2 + n+i] = -VR[j*n+i];
            A[(j+n) * n2 + n+i] = -VR[(j+n) * n + i];
        }


    zgeev("N", "V", &n2, A, &n2, W, VL, &n2, VR, &n2, work, &lwork, rwork, &info );

    if (info != 0)
    {
        printf ("error in zgeev in Krylov with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }



//  select modes with lamda_min<|lamda| < 1.0 + tol;
//  store the vector u to A matrix 
    num_modes = 0;
    for (i = 0; i < n2; i++) 
    {
        
        W[i] = (1.0/W[i] + shift);
        if(abs(W[i]) <1.0 - tol && abs(W[i]) > lamda_min)  // right going evanescent mode
        {
            
            lamda[num_modes] = W[i];
            zcopy(&n, &VR[i * n2], &ione, &A[num_modes*n], &ione);
            num_modes++;

        } 
        else if (abs(W[i]) >= 1.0-tol && abs(W[i]) < 1.0 + tol)  // propagating mode, calculate velocity
        {

            zgemv ("N", &n, &n, &W[i], H01, &n, &VR[i*n2], &ione, &zero, &VR[i*n2 + n], &ione);

            alpha = -zdotc(&n, &VR[i*n2], &ione, &VR[i*n2+n], &ione);
            if(alpha.imag() < 0.0)
            {
                lamda[num_modes] = W[i];
                zcopy(&n, &VR[i * n2], &ione, &A[num_modes*n], &ione);
                num_modes++;
            } 
        }

    }

//for (i = 0; i < n2; i++) 
//        printf("\n  %32.16f  lamda   %f  %f ", abs(W[i]), W[i]);
    //  now A store eigenvector n * num_modes  u(+) from literatures.
    //  find the dual vector u~ 
    // VR store overlap matrix S(num_modes, num_modes) 
    alpha = 1.0;
    zgemm("C", "N", &num_modes, &num_modes, &n, (double *)&alpha, (double *)A, &n, (double *)A, &n, (double *)&zero, (double *)VR, &num_modes);

    //  &A[num_modes *n] store transpose( u(+) )
    for(i = 0; i < num_modes; i++)
        for(j = 0; j < n; j++)
            A[num_modes * n + j * num_modes + i] = std::conj(A[i * n + j]);

    zgesv(&num_modes, &n, (double *)VR, &num_modes, ipiv, (double *)&A[num_modes*n], &num_modes, &info );
    if (info != 0)
    {
        printf ("error in zgesv in dual vector with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }

    for(i = 0; i < num_modes; i++)
        for(j = 0; j < n; j++)
        {
            A[ i * n + j] = lamda[i] * A[ i *n + j];
        }


    //  VR store F(+)
    alpha = 1.0;
    zgemm("N", "N", &n, &n, &num_modes, (double *)&alpha, (double *)A, &n, (double *)&A[num_modes*n], &num_modes, (double *)&zero, (double *)VR, &n);


    alpha = -1.0;
    zgemm("N", "N", &n, &n, &n, (double *)&alpha, (double *)H01, &n, (double *)VR, &n, (double *)&zero, (double *)sigma, &n);


    delete [] A;
    delete [] VR;
    delete [] W;
    delete [] lamda;
    delete [] work;
    delete [] rwork;
    delete [] ipiv;


}

