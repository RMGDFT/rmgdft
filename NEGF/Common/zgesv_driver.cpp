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
#include "blacs.h"


//#include "transition.h"



// Ax = B

void zgesv_driver (std::complex<double> *A, int *desca,  std::complex<double> *B, int *descb )
{



    int d_ipiv, *ipiv, info, ione =1, izero = 0;
    int nprow, npcol, myrow, mycol;
    int ictxt = desca[1];
    int nn = desca[2], mb= desca[4];
    int nhrs = descb[3];
    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    RmgTimer *RT5 = new RmgTimer("Ax=B");


    if(nprow*npcol <1) 
    {
        printf ("error in zgesv_driver nprow= %d npcol=%d \n", nprow, npcol);
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


        ipiv = (int *) malloc(d_ipiv * sizeof(int));
        DeviceSynchronize();
        magma_zgesv_gpu (nn, nhrs, (magmaDoubleComplex *)A, nn, ipiv, (magmaDoubleComplex *)B, nn, &info);

        if (info != 0)
        {
            printf ("error in magma_zgesv with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }
        free(ipiv);
    #else


        DeviceSynchronize();
        cusolverStatus_t cu_status;
        int Lwork;
        int *devIpiv, *devInfo;
        cuDoubleComplex *Workspace; 
        cudaError_t cuerr = cudaMalloc((void **)&devIpiv, sizeof(int) *nn);
        cuerr = cudaMalloc((void **)&devInfo, sizeof(int) );

        cu_status = cusolverDnZgetrf_bufferSize(ct.cusolver_handle, nn, nn, (cuDoubleComplex *)A, nn, &Lwork);
        if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__,"cusolverDnZgetrf_bufferSize failed.");
        cuerr = cudaMalloc((void **) &Workspace, sizeof(cuDoubleComplex) *Lwork);
        cu_status = cusolverDnZgetrf(ct.cusolver_handle, nn, nn, (cuDoubleComplex *)A, nn, Workspace, devIpiv, devInfo );
        if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__,"cusolverDnZgetrf failed.");
        info = 0;
        if (info != 0)
        {
            printf ("error in cusolverDnZgetrf with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }


        DeviceSynchronize();
        if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__,"cusolverDnZgetrf failed.");
        
        cublasOperation_t trans =CUBLAS_OP_N;
        cu_status = cusolverDnZgetrs(ct.cusolver_handle, trans, nn, nhrs, (const cuDoubleComplex *)A, nn, devIpiv, (cuDoubleComplex *)B, nn, devInfo ); 
        if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__,"cusolverDnZgetrs failed.");


        info = 0;
        if (info != 0)
        {
            printf ("error in cusolverDnZgetrs with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }
        //if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnZgetrs failed.");
        DeviceSynchronize();
        cudaFree(devIpiv);
        cudaFree(devInfo);
        cudaFree(Workspace);

    #endif

    
#else
    //  use scalapack if nprow * npcol > 1
    if(nprow*npcol > 1)  
    {

        d_ipiv = mb + numroc_( &nn, &mb, &myrow, &izero, &nprow);

        ipiv = (int *) malloc(d_ipiv * sizeof(int));
        pzgesv (&nn, &nhrs, A, &ione, &ione, desca, ipiv, B, &ione, &ione, descb, &info); 
        if (info != 0)
        {
            printf ("error in pzgesv with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }

        free(ipiv);

    }
    else
    {

        d_ipiv = nn;
        ipiv = (int *) malloc(d_ipiv * sizeof(int));

        zgesv(&nn, &nhrs, (double *)A, &nn, ipiv, (double *)B, &nn, &info );
        if (info != 0)
        {
            printf ("error in zgesv with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }

        free(ipiv);
    }
#endif

    delete RT5;

}

