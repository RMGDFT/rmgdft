#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

#include "RmgGemm.h"


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#include "blas.h"
#include "blacs.h"
#include "Scalapack.h"

#include "typedefs.h"
#include "rmg_control.h"
#include "transition.h"

#if HIP_ENABLED
#include <hip/hip_runtime.h>
#endif

#include "Gpufuncs.h"

void my_sync_device()
{
#if CUDA_ENABLED
    DeviceSynchronize();
#endif
#if HIP_ENABLED
    hipDeviceSynchronize();
#endif
}

void zcopy_driver (int n, std::complex<double> *A, int ia, std::complex<double> *B, int ib) 
{

#if CUDA_ENABLED 
    cublasZcopy (ct.gpublas_handle, n, (cuDoubleComplex *)A, ia, (cuDoubleComplex *)B, ib);
#elif HIP_ENABLED
    hipblasZcopy (ct.gpublas_handle, n, (hipblasDoubleComplex *)A, ia, (hipblasDoubleComplex *)B, ib);
#else
    zcopy (&n, A, &ia, B, &ib);
#endif
}


void zaxpy_driver (int n, std::complex<double> alpha, std::complex<double> *A, int ia, std::complex<double> *B, int ib) 
{

#if CUDA_ENABLED 
    cublasZaxpy (ct.gpublas_handle, n, (cuDoubleComplex *)&alpha, (cuDoubleComplex *)A, ia, (cuDoubleComplex *)B, ib);
#elif HIP_ENABLED
    hipblasZaxpy (ct.gpublas_handle, n, (hipblasDoubleComplex *)&alpha, (hipblasDoubleComplex *)A, ia, (hipblasDoubleComplex *)B, ib);
#else
    zaxpy (&n, &alpha, A, &ia, B, &ib);
#endif
}

void dzasum_driver(int n, std::complex<double> *A, int ia, double *sum)
{
#if CUDA_ENABLED 
    cublasDzasum (ct.gpublas_handle, n, (cuDoubleComplex *)A, ia, sum);
#elif HIP_ENABLED
    hipblasDzasum (ct.gpublas_handle, n, (hipblasDoubleComplex *)A, ia, sum);
#else
    *sum = dzasum(&n, (double *)A, &ia);
#endif
}


void dcopy_driver (int n, double *A, int ia, double *B, int ib) 
{

    if(!ct.tddft_gpu)
    {
        dcopy (&n, A, &ia, B, &ib);
    }
    else
    {
#if CUDA_ENABLED || HIP_ENABLED
        gpublasDcopy (ct.gpublas_handle, n, A, ia, B, ib);
#endif
    }
}

void scopy_driver (int n, float *A, int ia, float *B, int ib) 
{

#if CUDA_ENABLED || HIP_ENABLED
    gpublasScopy (ct.gpublas_handle, n, A, ia, B, ib);
#else
    scopy (&n, A, &ia, B, &ib);
#endif
}

void daxpy_driver (int n, double alpha, double *A, int ia, double *B, int ib) 
{

    if(!ct.tddft_gpu)
    {
        daxpy (&n, &alpha, A, &ia, B, &ib);
    }
    else
    {

#if CUDA_ENABLED || HIP_ENABLED
        gpublasDaxpy (ct.gpublas_handle, n, &alpha, A, ia, B, ib);
#endif
    }
}

void saxpy_driver (int n, float alpha, float *A, int ia, float *B, int ib) 
{

#if CUDA_ENABLED || HIP_ENABLED
    gpublasSaxpy (ct.gpublas_handle, n, &alpha, A, ia, B, ib);
#else
    saxpy (&n, &alpha, A, &ia, B, &ib);
#endif
}

void dscal_driver(int n, double beta, double *A, int ione)
{
    if(!ct.tddft_gpu)
    {
        dscal(&n, &beta, A, &ione);
    }
    else
    {
#if CUDA_ENABLED || HIP_ENABLED
        gpublasDscal (ct.gpublas_handle, n, &beta, A, ione);
#endif
    }

}

void dgemm_driver (char *transa, char *transb, int m, int n, int k, 
        double alpha, double *A, int ia, int ja, int *desca,
        double *B, int ib, int jb, int *descb, double beta, 
        double *C, int ic, int jc, int *descc)
{

    int nprow, npcol, myrow, mycol;
    int lda=desca[8], ldb=descb[8], ldc = descc[8];
    int ictxt = desca[1];

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);
    if(nprow*npcol <1) 
    {
        printf ("error in dgemmdriver nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }

    //  use scalapack if nprow * npcol > 1
    if(nprow*npcol > 1)  
    {
        pdgemm (transa, transb, &m, &n, &k, &alpha, A, &ia, &ja, desca,
                B, &ib, &jb, descb, &beta, C, &ic, &jc, descc);
    }
    else
    {
        RmgGemm (transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    }


}

void sgemm_driver (char *transa, char *transb, int m, int n, int k, 
        float alpha, float *A, int ia, int ja, int *desca,
        float *B, int ib, int jb, int *descb, float beta, 
        float *C, int ic, int jc, int *descc)
{

    int nprow, npcol, myrow, mycol;
    int lda=desca[8], ldb=descb[8], ldc = descc[8];
    int ictxt = desca[1];

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);
    if(nprow*npcol <1) 
    {
        printf ("error in sgemmdriver nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }

#if CUDA_ENABLED || HIP_ENABLED

    if(nprow*npcol != 1)
    {
        printf ("GPU ENABLED but nprow*npcol !=1  nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }

#endif

    //  use scalapack if nprow * npcol > 1
    if(nprow*npcol > 1)  
    {
        psgemm (transa, transb, &m, &n, &k, &alpha, A, &ia, &ja, desca,
                B, &ib, &jb, descb, &beta, C, &ic, &jc, descc);
    }
    else
    {
        RmgGemm (transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    }


}


void zgemm_driver (char *transa, char *transb, int m, int n, int k, 
        std::complex<double> alpha, std::complex<double> *A, int ia, int ja, int *desca,
        std::complex<double> *B, int ib, int jb, int *descb, std::complex<double> beta, 
        std::complex<double> *C, int ic, int jc, int *descc)
{

    int nprow, npcol, myrow, mycol;
    int lda=desca[8], ldb=descb[8], ldc = descc[8];
    int ictxt = desca[1];

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);
    if(nprow*npcol <1) 
    {
        printf ("error in zgemmdriver nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }

#if CUDA_ENABLED || HIP_ENABLED

    if(nprow*npcol != 1)
    {
        printf ("GPU ENALBED but nprow*npcol !=1  nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }

#endif

    //  use scalapack if nprow * npcol > 1
    if(nprow*npcol > 1)  
    {
        pzgemm (transa, transb, &m, &n, &k, &alpha, A, &ia, &ja, desca,
                B, &ib, &jb, descb, &beta, C, &ic, &jc, descc);
    }
    else
    {
        RmgGemm (transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    }


}


// Used for the special case when all MPI procs share a copy of the input matrices
// and the gemm call is repeated by all procs using the GPUs.
void mgpu_dgemm_driver (char *transa, char *transb, int m, int n, int k, 
        double alpha, double *A, int ia, int ja, int *desca,
        double *B, int ib, int jb, int *descb, double beta, 
        double *C, int ic, int jc, int *descc)
{

    int nprow, npcol, myrow, mycol;
    int lda=desca[8], ldb=descb[8], ldc = descc[8];
    int ictxt = desca[1];

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);
    if(nprow*npcol > 1)
    {
        dgemm_driver (transa, transb, m, n, k, 
                alpha, A, ia, ja, desca,
                B, ib, jb, descb, beta, 
                C, ic, jc, descc);
    }
    else
    {
#if CUDA_ENABLED || HIP_ENABLED

        if(m != n || m != k)
        {
            printf ("mgpu_dgemm requires m=n=k! \n");
            fflush (NULL);
            exit (0);
        }

#if USE_NCCL
        // gemm split over on node GPUs and recombined with nccl
        // local_comm contains the number of procs on this node which is also the number of GPUs
        int nprocs, my_rank;
        MPI_Comm_size(pct.local_comm, &nprocs);
        ncclCommUserRank(ct.nccl_local_comm, &my_rank);

        std::vector<int> start, stop, counts;
        start.resize(nprocs);
        stop.resize(nprocs);
        counts.resize(nprocs);

        int my_start = 0;
        int my_stop = 0;
        int my_step = 0;
        int incs = m / nprocs;
        if(m % nprocs) incs++;
        int ioffset = 0;
        for(int idx = 0;idx < nprocs;idx++)
        {
            start[idx] = ioffset;
            stop[idx] = start[idx] + incs;
            if(idx == (nprocs-1)) stop[idx] = m;
            counts[idx] = stop[idx] - start[idx];
            if(my_rank == idx)
            {
                my_start = start[idx];
                my_stop = stop[idx];
                my_step = counts[idx];
            }
            ioffset += counts[idx];
            start[idx] *= m;
            stop[idx] *= m;
            counts[idx] *= m;
        }
        RmgGemm(transa, transb, m, my_step, k, alpha, A, lda, &B[my_start*m], ldb, beta, &C[my_start*m], ldc);
        size_t sendcount = counts[0];
        ncclAllGather(&C[my_rank*sendcount], C, sendcount, ncclDouble, ct.nccl_local_comm, 0);
#else
        // Full gemm no nccl
        RmgGemm (transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
#endif

    }


}


