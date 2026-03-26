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

#if HIP_ENABLED
    #include <hipblas/hipblas.h>
    #include <rocsolver/rocsolver.h>
    #include <hip/hip_runtime_api.h> // for hip functions
    #include <hipsolver/hipsolver.h> // for all the hipsolver C interfaces and type declarations
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
#include "rmg_gemm.h"
#include "GpuAlloc.h"

#include "Gpufuncs.h"
#include "blas.h"
#include "Scalapack.h"
#include "GpuAlloc.h"
#include "RmgMatrix.h"
#include "transition.h"






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
        rmg::printlog ("error in matrix_inverse_driver nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }

#if CUDA_ENABLED

    if(nprow*npcol != 1)
    {
        rmg::printlog ("GPU ENALBED but nprow*npcol !=1  nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }
    
    std::complex<double> *gpu_temp, *Imatrix;
    size_t size = nn*nn*sizeof(std::complex<double>);
    Imatrix = (std::complex<double> *)RmgMallocHost(size);
    pmo_unitary_matrix(Imatrix, desca);
    gpuMalloc((void **)&gpu_temp, size);
    MemcpyHostDevice(size, Imatrix, gpu_temp);


    std::complex<double> *A = Hii;
    std::complex<double> *B = gpu_temp;
    
    rmg::sync_device();
    cusolverStatus_t cu_status;
    int Lwork;
    int *devIpiv, *devInfo;
    cuDoubleComplex *Workspace;
    gpuMalloc((void **)&devIpiv, sizeof(int) *nn);
    gpuMalloc((void **)&devInfo, sizeof(int) );

    cu_status = cusolverDnZgetrf_bufferSize(ct.cusolver_handle, nn, nn, (cuDoubleComplex *)A, nn, &Lwork);
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg::error("cusolverDnZgetrf_bufferSize failed.");
    gpuMalloc((void **) &Workspace, sizeof(cuDoubleComplex) *Lwork);
    cu_status = cusolverDnZgetrf(ct.cusolver_handle, nn, nn, (cuDoubleComplex *)A, nn, Workspace, devIpiv, devInfo );
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg::error("cusolverDnZgetrf failed.");
    info = 0;
    if (info != 0)
    {
        rmg::printlog ("error in cusolverDnZgetrf with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }


    rmg::sync_device();
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg::error("cusolverDnZgetrf failed.");

    cublasOperation_t trans =CUBLAS_OP_N;
    cu_status = cusolverDnZgetrs(ct.cusolver_handle, trans, nn, nn, (const cuDoubleComplex *)A, nn, devIpiv, (cuDoubleComplex *)B, nn, devInfo );
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg::error("cusolverDnZgetrs failed.");


    info = 0;
    if (info != 0)
    {
        rmg::printlog ("error in cusolverDnZgetrs with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }
    //if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg::error(" cusolverDnZgetrs failed.");
    rmg::sync_device();
    gpuFree(devIpiv);
    gpuFree(devInfo);
    gpuFree(Workspace);

    rmg::zcopy_driver (nn*nn, gpu_temp, ione, Hii, ione);
    RmgFreeHost(Imatrix);
    gpuFree(gpu_temp);

#elif HIP_ENABLED

    if(nprow*npcol != 1)
    {
        rmg::printlog ("GPU ENALBED but nprow*npcol !=1  nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }

    std::complex<double> *gpu_temp, *Imatrix;
    size_t size = nn*nn*sizeof(std::complex<double>);
    Imatrix = (std::complex<double> *)RmgMallocHost(size);
    pmo_unitary_matrix(Imatrix, desca);
    gpuMalloc((void **)&gpu_temp, size);
    MemcpyHostDevice(size, Imatrix, gpu_temp);


    ZgetrftrsDriver(nn, nn, Hii, gpu_temp);

    rmg::zcopy_driver (nn*nn, gpu_temp, ione, Hii, ione);
    RmgFreeHost(Imatrix);
    gpuFree(gpu_temp);

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
            rmg::printlog ("error in pzgetrf with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }
        pzgetri(&nn, Hii, &ione, &ione, desca, ipiv, work, &lwork, iwork, &liwork, &info);
        if (info != 0)
        {
            rmg::printlog ("error in pzgetri with INFO = %d \n", info);
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

        zgetrf(&nn, &nn, Hii, &nn, ipiv, &info);
        if (info != 0)
        {
            rmg::printlog ("error in zgetrf with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }
        zgetri(&nn, Hii, &nn, ipiv, work, &lwork, &info);
        if (info != 0)
        {
            rmg::printlog ("error in zgetri with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }

        free(ipiv);
        free(work);
    }
#endif

    delete RT5;

}

