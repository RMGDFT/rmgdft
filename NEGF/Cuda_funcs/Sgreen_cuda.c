/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
 *   Iterative construction of the transfer matrix, 
 *	as Lopez-Sancho & Rubio, J.Phys.F: Met. Phs., v.14, 1205 (1984) 
 *		and ibid. v.15, 851 (1985)
 *
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <complex.h>

#include "main.h"
#include "init_var_negf.h"
#include "LCR.h"
#include "pmo.h"

#define 	MAX_STEP 	100


void Sgreen_cuda (complex double * g, 
        complex double *ch0, complex double *ch01, complex double *ch10, int iprobe)
{

#if GPU_ENABLED
    cuDoubleComplex cumone, cuone,  cuzero, *Hiii;
    cublasOperation_t transT = CUBLAS_OP_T, transN = CUBLAS_OP_N;
    cuDoubleComplex *tau, *taut, *tsum, *tsumt, *t11, *t12, *s1, *s2;
    cuDoubleComplex *ch0_cu, *ch1_cu, *tot_cu, *tott_cu;
    complex double *Imatrix;
    int *ipiv;
    int n1, nmax, info, ione = 1, i, step;
    double converge1, converge2;

    nmax = lcr[iprobe].num_states;
    n1 = nmax * nmax;
    cuone.x = 1.0;
    cumone.x = -1.0;
    cuzero.x = 0.0;
    cuone.y = 0.0;
    cumone.y = 0.0;
    cuzero.y = 0.0;


    my_malloc_init( Imatrix, nmax * nmax, complex double );
    my_malloc_init( ipiv, nmax, int );


    for (i = 0; i < nmax * nmax; i++)
    {
        Imatrix[i] = 0.0;
    }

    for (i = 0; i < nmax; i++)
    {
        Imatrix[i * nmax + i] = 1.0;
    }


    tau = ct.gpu_Htri;
    taut = ct.gpu_Htri + 2 * nmax * nmax;
    tsum = ct.gpu_Gtri;
    tsumt = ct.gpu_Gtri + nmax * nmax;
    t11 = ct.gpu_Gtri + 2 * nmax * nmax;
    t12 = ct.gpu_Gtri + 3 * nmax * nmax;
    s1 = ct.gpu_Gtem;
    s2 = ct.gpu_Gtem + nmax * nmax;
    ch0_cu = ct.gpu_temp;
    ch1_cu = ct.gpu_Gtem + 2* nmax * nmax;
    tot_cu = ct.gpu_Hii;
    tott_cu = ct.gpu_Gii;

    cublasSetVector( n1, sizeof( complex double ), Imatrix, ione, ct.gpu_Imatrix, ione );
    cublasSetVector( n1, sizeof( complex double ), ch0, ione, ch0_cu, ione );

    cublasZscal (ct.cublas_handle, n1, &cuzero, t12, ione);
    cublasZaxpy (ct.cublas_handle, n1, &cumone, ch0_cu, ione, t12, ione);

    cublasZcopy (ct.cublas_handle, n1, ct.gpu_Imatrix, ione, t11, ione);



    magma_zgesv_gpu( nmax, nmax, t12, nmax, ipiv, t11, nmax, &info );


    if (info != 0)
    {
        printf ("Stransfer_cuda.c: error in PZGESV with INFO = %d in pe %d\n", info, pct.gridpe);
        fflush (NULL);
        MPI_Finalize ();
        exit (0);
    }
    fflush(NULL);

    /* initialize intermediate t-matrices  */

    cublasSetVector( n1, sizeof( complex double ), ch10, ione, ch1_cu, ione );
    cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, t11, nmax,
            ch1_cu, nmax, &cuzero, tau, nmax);

    cublasSetVector( n1, sizeof( complex double ), ch01, ione, ch1_cu, ione );
    cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, t11, nmax,
            ch1_cu, nmax, &cuzero, taut, nmax);

    cublasZcopy (ct.cublas_handle, n1, tau, ione, tot_cu, ione);
    cublasZcopy (ct.cublas_handle, n1, taut, ione, tott_cu, ione);
    cublasZcopy (ct.cublas_handle, n1, taut, ione, tsum, ione);
    cublasZcopy (ct.cublas_handle, n1, tau, ione, tsumt, ione);



    /*  iterative loop till convergence is achieved  */
    for (step = 0; step < MAX_STEP; step++)
    {


        cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, tau, nmax,
                taut, nmax, &cuzero, t11, nmax);
        cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, taut, nmax,
                tau, nmax, &cuzero, t12, nmax);

        cublasZcopy (ct.cublas_handle, n1, ct.gpu_Imatrix, ione, s1, ione);
        cublasZcopy (ct.cublas_handle, n1, ct.gpu_Imatrix, ione, s2, ione);

        cublasZaxpy (ct.cublas_handle, n1, &cumone, t11, ione, s1, ione);
        cublasZaxpy (ct.cublas_handle, n1, &cumone, t12, ione, s1, ione);


        magma_zgesv_gpu( nmax, nmax, s1, nmax, ipiv, s2, nmax, &info );
        if (info != 0)
        {
            printf ("Stransfer_cuda.c: error in ZGESV with INFO = %d %d\n", info, step);
            fflush (NULL);
            MPI_Finalize ();
            exit (0);
        }

        cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, tau, nmax,
                tau, nmax, &cuzero, t11, nmax);
        cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, taut, nmax,
                taut, nmax, &cuzero, t12, nmax);
        cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, s2, nmax,
                t11, nmax, &cuzero, &tau[n1], nmax);
        cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, s2, nmax,
                t12, nmax, &cuzero, &taut[n1], nmax);
        cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, tsum, nmax,
                &tau[n1], nmax, &cuzero, t11, nmax);
        cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, tsum, nmax,
                &taut[n1], nmax, &cuzero, s1, nmax);

        cublasZcopy (ct.cublas_handle, n1, t11, ione, s2, ione);

        cublasZaxpy (ct.cublas_handle, n1, &cuone, tot_cu, ione, s2, ione);

        cublasZcopy (ct.cublas_handle, n1, s2, ione, tot_cu, ione);
        cublasZcopy (ct.cublas_handle, n1, s1, ione, tsum, ione);

        cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, tsumt, nmax,
                &taut[n1], nmax, &cuzero, t11, nmax);
        cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, tsumt, nmax,
                &tau[n1], nmax, &cuzero, s1, nmax);

        cublasZcopy (ct.cublas_handle, n1, t11, ione, s2, ione);

        cublasZaxpy (ct.cublas_handle, n1, &cuone, tott_cu, ione, s2, ione);


        cublasZcopy (ct.cublas_handle, n1, s2, ione, tott_cu, ione);
        cublasZcopy (ct.cublas_handle, n1, s1, ione, tsumt, ione);
        cublasZcopy (ct.cublas_handle, n1, &tau[n1], ione, tau, ione);
        cublasZcopy (ct.cublas_handle, n1, &taut[n1], ione, taut, ione);


        cublasDzasum(ct.cublas_handle, n1, &tau[n1], ione, &converge1);
        cublasDzasum(ct.cublas_handle, n1, &taut[n1], ione, &converge2);


        if (converge1 < 1.0E-7 && converge2 < 1.0E-7)
            break;
    }

    if (converge1 > 1.0E-7 || converge2 > 1.0E-7)
    {
        printf ("bad t-matrix convergence\n");
        fflush (NULL);
        MPI_Finalize ();
        exit (0);
    }
//    cublasGetVector(n1, sizeof( complex double ), tot_cu, ione, tot, ione );
//    cublasGetVector(n1, sizeof( complex double ), tott_cu, ione, tott, ione );

        cublasZgemm (ct.cublas_handle, transN, transN, nmax, nmax, nmax, &cuone, ch1_cu, nmax,
                tot_cu, nmax, &cuone, ch0_cu, nmax);

        magma_zgesv_gpu( nmax, nmax, ch0_cu, nmax, ipiv, ct.gpu_Imatrix, nmax, &info );
    if (info != 0)
    {
        dprintf ("Sgreen_cuda.c: error in last ZGESV  with INFO = %d \n", info);
        fflush (NULL);
        MPI_Finalize ();
        exit (0);
    }

    cublasGetVector(n1, sizeof( complex double ), ct.gpu_Imatrix, ione, g, ione );

    my_free(ipiv);
    my_free( Imatrix );
#endif

}


