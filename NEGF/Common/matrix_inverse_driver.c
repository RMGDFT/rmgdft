/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "make_conf.h"

#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "LCR.h"
#include "init_var.h"

#if GPU_ENABLED
#include <magma.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#include "Scalapack.h"
#include "blas.h"
#include "RmgTimer.h"



void matrix_inverse_driver (complex double *Hii, int *desca )
{



    int d_ipiv, *ipiv, info, ione =1, izero = 0;
    int nprow, npcol, myrow, mycol;
    int lwork, liwork, *iwork;
    complex double *work;
    int ictxt = desca[1];
    int nn = desca[2], mb= desca[4], nb=desca[5];
    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    void *RT5 = BeginRmgTimer("matrix inverse");


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


    magma_zgetrf_gpu(nn, nn, (magmaDoubleComplex *)Hii, nn, ipiv, &info);
    if (info != 0)
    {
        printf ("error in magma_zgetrf with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }
    magma_zgetri_gpu(nn, (magmaDoubleComplex *)Hii, nn, ipiv, (magmaDoubleComplex *)ct.gpu_temp, lwork, &info);
    if (info != 0)
    {
        printf ("error in magma_zgetri with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }
    free(ipiv);
#else
    //  use scalapack if nprow * npcol > 1
    if(nprow*npcol > 1)  
    {

        d_ipiv = mb + numroc_( &nn, &mb, &myrow, &izero, &nprow);
        double work_tem;
        int ltem = -1;
        pzgetri_(&nn, (double *)Hii, &ione, &ione, desca, ipiv, &work_tem, &ltem,  &liwork, &ltem, &info);

        lwork = (int)work_tem +1;

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
            printf ("error in zgetrf with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }
        zgetri(&nn, (double *)Hii, &nn, ipiv, (double *)work, &lwork, &info);
        if (info != 0)
        {
            printf ("error in zgetri with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }

        free(ipiv);
        free(work);
    }
#endif

    EndRmgTimer(RT5);

}

