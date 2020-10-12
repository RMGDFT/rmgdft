#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if CUDA_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
    #include <cusolverDn.h>
    #if MAGMA_LIBS
        #include <magma.h>
    #endif
#endif


#include <complex>

#include "Scalapack.h"
#include "blas.h"
#include "RmgTimer.h"

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "Subdiag.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "Gpufuncs.h"
#include "blas.h"
#include "Scalapack.h"
#include "GpuAlloc.h"





void matrix_inverse_driver (std::complex<double> *Hii, int *desca )
{



    int d_ipiv, *ipiv=NULL, info, ione =1, izero = 0;
    int nprow, npcol, myrow, mycol;
    int lwork, liwork, *iwork;
    std::complex<double> *work;
    int ictxt = desca[1];
    int nn = desca[2], mb= desca[4];
    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    RmgTimer *RT5 = new RmgTimer("matrix inverse");


    if(nprow*npcol <1) 
    {
        printf ("error in matrix_inverse_driver nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }
#if CUDA_ENABLED

    if(nprow*npcol != 1)
    {
        printf ("GPU ENALBED but nprow*npcol !=1  nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }
    
#if MAGMA_LIBS
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
    
    std::complex<double> *gpu_temp = (std::complex<double> *)GpuMallocManaged(nn*nn*sizeof(std::complex<double>));
    pmo_unitary_matrix(gpu_temp, desca);
    zgesv_driver (Hii, desca, gpu_temp, desca);
    zcopy_driver (nn*nn, gpu_temp, ione, Hii, ione);
    GpuFreeManaged(gpu_temp);

#endif

#else
    //  use scalapack if nprow * npcol > 1
    if(nprow*npcol > 1)  
    {

        d_ipiv = mb + numroc_( &nn, &mb, &myrow, &izero, &nprow);
        std::complex<double> work_tem;
        int ltem = -1;
        pzgetri(&nn, Hii, &ione, &ione, desca, ipiv, &work_tem, &ltem,  &liwork, &ltem, &info);

        lwork = (int)std::real(work_tem) +1;

        ipiv = (int *) malloc(d_ipiv * sizeof(int));
        iwork = (int *) malloc(liwork * sizeof(int));
        work = (std::complex<double> *)malloc(lwork * sizeof(std::complex<double>));

        pzgetrf(&nn, &nn, Hii, &ione, &ione, desca, ipiv, &info);
        if (info != 0)
        {
            printf ("error in pzgetrf with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }
        pzgetri(&nn, Hii, &ione, &ione, desca, ipiv, work, &lwork, iwork, &liwork, &info);
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
        work = (std::complex<double> *)malloc(lwork * sizeof(std::complex<double>));

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

    delete RT5;

}

