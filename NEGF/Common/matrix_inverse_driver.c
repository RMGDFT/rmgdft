/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "make_conf.h"

#if GPU_ENABLED
//#include <magma.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#include "my_scalapack.h"
#include "blas.h"



void matrix_inverse_driver (complex double *Hii, int *desca )
{



    int d_ipiv, *ipiv, info, ione =1, izero = 0;
    int nprow, npcol, myrow, mycol;
    int lwork, liwork, *iwork;
    complex double *work;
    int ictxt = desca[1];
    int nn = desca[2], mb= desca[4], nb=desca[5];
    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);


    if(nprow*npcol <1) 
    {
        printf ("error in matrix_inverse_driver nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }
#if GPU_ENABLED

    if(nprow*npcol != 1)
    {
        printf ("GPU ENALBED but nprow*npcol !=1  nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }
    
    d_ipiv = nn;
    lwork = magma_get_zgetri_nb(nn);
    lwork = lwork * nn;
    //lwork = nn * nn;

    cudaError_t cuerr;


    ipiv = (int *) malloc(d_ipiv * sizeof(int));
    cuerr = cudaMalloc((void **)&work , lwork * sizeof(complex double) );
    if(cuerr != cudaSuccess)
    {
        printf ("cuda allocation faild in matrix_inverse_driver lwork = %d \n", lwork);
        fflush (NULL);
        exit (0);
    }


    magma_zgetrf_gpu(nn, nn, Hii, nn, ipiv, &info);
    if (info != 0)
    {
        printf ("error in magma_zgetrf with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }
    magma_zgetri_gpu(nn, Hii, nn, ipiv, work, lwork, &info);
    if (info != 0)
    {
        printf ("error in magma_zgetri with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }
    free(ipiv);
    cudaFree(work);
#else
    //  use scalapack if nprow * npcol > 1
    if(nprow*npcol > 1)  
    {

        d_ipiv = mb + numroc_( &nn, &mb, &myrow, &izero, &nprow);
        lwork = d_ipiv * nb;
        liwork = d_ipiv + mb + ceil(nn/mb);
        ipiv = (int *) malloc(d_ipiv * sizeof(int));
        iwork = (int *) malloc(liwork * sizeof(int));
        work = (complex double *)malloc(lwork * sizeof(complex double));

        pzgetrf_(&nn, &nn, (double *)Hii, &ione, &ione, desca, ipiv, &info);
        if (info != 0)
        {
            printf ("error in pzgetrf with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }
        pzgetri_(&nn, (double *)Hii, &ione, &ione, desca, ipiv, (double *)work, &lwork, iwork, &liwork, &info);
        if (info != 0)
        {
            printf ("error in pzgetri with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }

        free(ipiv);
        free(iwork);
        free(work);

    }
    else
    {


        d_ipiv = nn;
        lwork = nn * nn;
        ipiv = (int *) malloc(d_ipiv * sizeof(int));
        work = (complex double *)malloc(lwork * sizeof(complex double));

        zgetrf(&nn, &nn, (double *)Hii, &nn, ipiv, &info);
        if (info != 0)
        {
            printf ("error in pzgetrf with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }
        zgetri(&nn, (double *)Hii, &nn, ipiv, (double *)work, &lwork, &info);
        if (info != 0)
        {
            printf ("error in pzgetri with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }

        free(ipiv);
        free(work);
    }
#endif


}

