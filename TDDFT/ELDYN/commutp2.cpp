  
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "blas.h"
#include "blacs.h"
#include "blas_driver.h"
#include "Scalapack.h"
#include "GpuAlloc.h"
#include "RmgTimer.h"
#include "main.h"

#include <iostream>
//#include "lapacke.h"
#if CUDA_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include "GpuAlloc.h"
#endif

#if HIP_ENABLED
    #include <hipblas/hipblas.h>
    #include "GpuAlloc.h"
#endif


#include "rmg_gemm.h"

#include "Gpufuncs.h"
#include "rmg_error.h"
#include "rmg_hvector.h"
#include "rmg_dvector.h"



void print_matrix(double *, int ) ;
void print_matrix2(double *, int ) ;

/////////////////////////////////////////////////////////////////////////
void TiledM_transpose(double *A, int numst)
{
#if CUDA_ENABLED || HIP_ENABLED

    rmg::error("\n not work yet\n");

#else
   rmg::hvector<double> A_glob(numst*numst);
   TiledM_to_glob(A_glob.data(), A, numst, pct.local_comm);
   int numst_pe = numst/pct.local_comm_npes;
   for(int i = 0; i < numst_pe; i++)
   {
       for(int j = 0; j < numst; j++)
       {
           A[i * numst + j] = A_glob[j * numst + i + pct.local_rank *numst_pe];
       }
   }

#endif
}

void transpose( double *A,  double *B, int *desca) 
{
    //printf(" transpose  %d \n", Nbasis) ;
    int Nbasis = desca[3];
    int ictxt = desca[1], nprow, npcol, myrow, mycol;

    rmg::sync_device();

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);
    if(nprow*npcol <1) 
    {
        printf("\n nprow npcol not right in transpose %d %d\n", nprow, npcol);
        fflush(NULL);
        exit(0);
    }
    else if(nprow * npcol == 1)
    {
#if CUDA_ENABLED || HIP_ENABLED
        double alpha = 1.0, beta = 0.0;
        rmg::error(gpublasDgeam(ct.gpublas_handle, GPUBLAS_OP_T, GPUBLAS_OP_N, Nbasis, Nbasis, &alpha,
                    A, Nbasis, &beta, B, Nbasis, B, Nbasis));

#else
        int ij=0  ;
        for (int ix=0; ix< Nbasis;ix++) {
            for (int iy=0; iy<Nbasis;  iy++) {
                B[iy*Nbasis +ix] =A[ij++]  ;
            }   
        }
#endif
    }
    else
    {
        double zero = 0.0, one = 1.0;
        int ione = 1;
        pdtran(&Nbasis, &Nbasis, &one, A, &ione, &ione, desca, &zero, B, &ione, &ione, desca);
    }
    //printf(" *** A :     \n"); print_matrix(A,Nbasis);
    //printf(" *** B= A^T :\n"); print_matrix(B,Nbasis);
}


/////////////////////////////////////////////////////////////////////////

void  commstr(double *A, double *B, double *C,double *p_alpha,int *desca, int Mdim, int Ndim, double *W)
{
    /*\
      ! commutator  of two symmetric matrices:
      !   [S1,S2] =S1*S2-*S2*S1= S1*S2 -transpose(S1*S2)
      \*/
    double  alpha  = *(p_alpha)  ;
    double  beta   = 0.0e0       ;
    int     Nbsq   = Mdim * Ndim ;
    int Mglob = desca[3];
    int ione = 1;

    char trN='N' ;  // transpose ='not
                    // call dgemm('N','N',Nb,Nb,Nb,alpha,A,Nb,B,Nb,beta,C,Nb)
    rmg::sync_device();
    RmgTimer *RT = new RmgTimer("2-TDDFT: dgemm");
    rmg::mgpu_dgemm_driver (&trN, &trN, Mglob, Mglob, Mglob, alpha, A, ione, ione, desca,
            B, ione, ione, desca, beta, C, ione, ione, desca);
    rmg::sync_device();
    delete RT;
    RT = new RmgTimer("2-TDDFT: transpose");
    transpose(C,W,desca) ;
    rmg::sync_device();
    delete RT;

    //SUBROUTINE DAXPY    ( n, alpha, x, incx, y, incy ) :: //  y  <--  alpha*x + y
    double minus_one = -1.0e0 ;

    rmg::daxpy_driver(Nbsq,  minus_one , W, ione, C, ione) ;
}

/////////////////////////////////////////////////////////////////////////
void  commatr(double *A, double *B, double *C,double *p_alpha,int *desca, int Mdim, int Ndim, double *W)
{
    /*\
      ! commutator of two anti-symmetric matrices:
      !  [S,A] = S*A-A*S =  S*A +transpose(S*A)
      \*/
    double  alpha  = *(p_alpha)  ;
    double  beta   = 0.0e0       ;
    int     Nbsq   = Mdim * Ndim ;

    char trN='N' ;  // transpose ='not
                    // call dgemm('N','N',Nb,Nb,Nb,alpha,A,Nb,B,Nb,beta,C,Nb)
    int ione = 1, Mglob = desca[3];
    rmg::sync_device();
    RmgTimer *RT = new RmgTimer("2-TDDFT: dgemm");
    rmg::mgpu_dgemm_driver (&trN, &trN, Mglob, Mglob, Mglob, alpha, A, ione, ione, desca,
            B, ione, ione, desca, beta, C, ione, ione, desca);
    delete RT;
    rmg::sync_device();
    RT = new RmgTimer("2-TDDFT: transpose");
    transpose(C,W,desca) ;
    rmg::sync_device();
    delete RT;

    double one = 1.0e0 ;
    rmg::daxpy_driver(Nbsq,  one , W, ione, C, ione) ;
}


//////////////////////////////////////////////////////////////////////////////////////

void  tstconv(double *C,int *p_M, double *p_thrs,int *p_ierr, double *p_err, bool *p_tconv, MPI_Comm comm) 
{

    int     M     = *(p_M)    ;  //  [in]  :  total  size of matrix (2*Nbasis*Nbasis)
    double  thrs  = *(p_thrs) ;  //  [in]  :  convergence threshold

    double  err               ;  //  [out] :   error= abs of max element in the matrix
    int    ierr   =   0       ;  //  [out] :   location of err in matrix/vector
    bool   tconv  =  false    ;  //  [out] :   if converged ?  true or false?


    rmg::sync_device();
    if(!ct.tddft_gpu)
    { 
        err = abs(C[0]); 
        for (int i=0; i <M ;i++) {
            double err_tmp = abs(C[i]) ; 
            if (err_tmp > err) {
                err  = err_tmp ;
                ierr = i       ;
            }
        }
    }
    else
    {
#if CUDA_ENABLED || HIP_ENABLED
        int idx;
#if HIP_ENABLED
        hipblasIdamax(ct.gpublas_handle, M, C, 1, &idx);
#endif
#if CUDA_ENABLED
        cublasIdamax(ct.gpublas_handle, M, C, 1, &idx);
#endif
        idx -=1;    
        // hipblasIdamax return the index in fortran way, starting from 1
        gpuMemcpy(&err, &C[idx], sizeof(double), gpuMemcpyDeviceToHost);
        err = abs(err);
        ierr = idx;
#endif
    }

    rmg::sync_device();
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, comm);

    if (err < thrs)  tconv = true ;
    /*-- return values **/
    *(p_err )  =  err  ;
    *(p_ierr)  =  ierr ;
    *(p_tconv) =  tconv ;
}


//////////////////////////////////////////////////////////////////////////////////////
void commutp(double *P0, double *P1, double *Om, int *desca, int Mdim, int Ndim,
        double *p_thrs, int *p_maxiter, int *p_niter, double *p_errmax, int *p_iprint, MPI_Comm comm) 
{
    /*\
      !
      ! P = P0 +(dt/ih) [F,P]+ (dt/ih)^2/2! [F[F,P] + (dt/ih)^3/3! [F,F[F,P]] +... 
      !
      ! P(k)  <-- P(k-1) + 1/(ihk) *[Om,dP(k-1)] = 
      !         = P(k-1) +1/(hk)  *[Om,A ] -i/(hk) *[Om,S ]  
      !   
      !   where dP = S+i*A 
      !  
      !  Om = (Fav*dt)/I   from Magnus
      !
      \*/
    double   thrs    = *(p_thrs)       ;
    double   maxiter = *(p_maxiter)    ;
    int      iprint  = *(p_iprint)     ;

    int      Nsq = Mdim * Ndim         ;
    int      Nsq2 = 2*Nsq              ;
    int      iter   = 1                ;
    bool     tConv  = false            ;
    double   errmax = 0.0e0            ;
    double   err    = 0.0e0            ;
    int      ierr   = 0                ;
    int      ione   = 1                ;
    double   rone   = 1.0e0            ;
    double   mone   = -1.0e0           ;
    double   zero = 0.0;

    if(ct.tddft_tiledMM)
    {
        int Mglob = std::max(Mdim, Ndim);
        double *Om_glob, *C, *dP;
        MallocHostOrDevice((void **)&Om_glob,  sizeof(double) * Mglob * Mglob);
        MallocHostOrDevice((void **)&C,  sizeof(double) * Nsq2);
        MallocHostOrDevice((void **)&dP,  sizeof(double) * Nsq2);

        TiledM_to_glob(Om_glob, Om, Mglob, pct.local_comm);

        rmg::dcopy_driver(Nsq2,  P0, ione, P1   ,ione) ;   // P1 =P0:  saves P0 into P1 (for updates) 
        rmg::dcopy_driver(Nsq2,  P0, ione, dP   ,ione) ;   // dP =P0:  saves P0 into dP (for commutator )
                                                           //
        while ( iter  <= maxiter && tConv ==  false ) {
            double alpha     =  1.0e0 /iter ;
            double neg_alpha = -alpha       ;
            // dP first Nsq elements are real part and second Nsq elements are imaginary part.
            //
            rmg::gemm("N", "N", Mglob, Ndim, Mglob, alpha, Om_glob, Mglob, &dP[Nsq], Mglob, zero, C, Mglob);
            rmg::gemm("N", "N", Mglob, Ndim, Mglob, neg_alpha, Om_glob, Mglob, &dP[0], Mglob, zero, &C[Nsq], Mglob);
            rmg::dcopy_driver(Nsq2, C, ione, dP   ,ione)             ;  

            TiledM_transpose(C, Mglob);
            rmg::daxpy_driver(Nsq,  rone , C, ione, dP, ione) ;

            TiledM_transpose(&C[Nsq], Mglob);
            rmg::daxpy_driver(Nsq,  mone , &C[Nsq], ione, &dP[Nsq], ione) ;

            rmg::daxpy_driver(Nsq2, rone,  dP, ione, P1, ione)       ;  // P1 =P1 +dP
            tstconv(dP, &Nsq2, &thrs,&ierr,&err,&tConv, comm)  ;  // tstconv(dP,2*Nsq,N,thrs,ierr,err,tconv)
            if (iprint>0) rmg::printlog("ConvergTest: Niter  %d  errmax = %10.5e \n",  iter,err) ;
            if (abs(err) >  errmax)  errmax= abs(err)  ;

            iter ++ ;

        }


        FreeHostOrDevice(dP);
        FreeHostOrDevice(C);
        FreeHostOrDevice(Om_glob);
        return;
    }


    // Needed for padding in C_dev
    if(!ct.tddft_gpu)
    {
        /* ----  P1=P0 ,  dP=P0   ---- */
        double      *C   = (double *)RmgMallocHost(Nsq2 * sizeof(double));
        double      *dP  = (double *)RmgMallocHost(Nsq2 * sizeof(double));
        rmg::dcopy_driver(Nsq2,  P0, ione, P1   ,ione) ;   // P1 =P0:  saves P0 into P1 (for updates) 
        rmg::dcopy_driver(Nsq2,  P0, ione, dP   ,ione) ;   // dP =P0:  saves P0 into dP (for commutator )
        double *W   = (double *)RmgMallocHost(Mdim * Ndim * sizeof(double));
        while ( iter  <= maxiter && tConv ==  false ) {
            double alpha     =  1.0e0 /iter ;
            double neg_alpha = -alpha       ;

            commatr(Om, &dP[Nsq], &C[0]  , &alpha,    desca, Mdim, Ndim, W)    ;  // real contrib correction  
            commstr(Om, &dP[0]  , &C[Nsq], &neg_alpha,desca, Mdim, Ndim, W)    ;  // imag contrib correction 

            rmg::dcopy_driver(Nsq2, C, ione, dP   ,ione)             ;  // dP=C     
            rmg::daxpy_driver(Nsq2, rone,  dP, ione, P1, ione)       ;  // P1 =P1 +dP
            tstconv(dP, &Nsq2, &thrs,&ierr,&err,&tConv, comm)  ;  // tstconv(dP,2*Nsq,N,thrs,ierr,err,tconv)
            if (iprint>0) rmg::printlog("ConvergTest: Niter  %d  errmax = %10.5e \n",  iter,err) ;
            if (abs(err) >  errmax)  errmax= abs(err)  ;

            iter ++ ;
        }
        RmgFreeHost(W);
        RmgFreeHost(C);
        RmgFreeHost(dP);
    }
    else
    {
#if CUDA_ENABLED || HIP_ENABLED
        int      Nsq_alloc                 ;
        Nsq_alloc = (Mdim + pct.procs_per_host)*Ndim;
        double *Om_dev, *P0_dev, *P1_dev, *W_dev, *C_dev, *dP_dev;
        gpuMalloc((void **)&Om_dev, Nsq * sizeof(double));
        gpuMalloc((void **)&P0_dev, Nsq2 * sizeof(double));
        gpuMalloc((void **)&P1_dev, Nsq2 * sizeof(double));
        gpuMalloc((void **)&W_dev, Nsq * sizeof(double));
        //    gpuMalloc((void **)&C_dev, Nsq2_alloc * sizeof(double));
        double *C0_dev, *C1_dev;
        gpuMalloc((void **)&C0_dev, Nsq_alloc * sizeof(double));
        gpuMalloc((void **)&C1_dev, Nsq_alloc * sizeof(double));
        gpuMalloc((void **)&dP_dev, Nsq2 * sizeof(double));
        rmg::dcopy_driver(Nsq,   Om, ione, Om_dev ,ione) ;      // GPU buffer for Om
        /* ----  P1=P0 ,  dP=P0   ---- */
        rmg::dcopy_driver(Nsq2,  P0, ione, P1_dev ,ione) ;   // P1 =P0:  saves P0 into P1 (for updates) 
        rmg::dcopy_driver(Nsq2,  P0, ione, dP_dev, ione) ;   // dP =P0:  saves P0 into dP (for commutator )
        while ( iter  <= maxiter && tConv ==  false ) {
            double alpha     =  1.0e0 /iter ;
            double neg_alpha = -alpha       ;

            //commatr(Om_dev, &dP_dev[Nsq], &C_dev[0]  , &alpha,    desca, Mdim, Ndim, W_dev)    ;  // real contrib correction  
            commatr(Om_dev, &dP_dev[Nsq], C0_dev  , &alpha,    desca, Mdim, Ndim, W_dev)    ;  // real contrib correction  
                                                                                               //        commstr(Om_dev, &dP_dev[0]  , &C_dev[Nsq], &neg_alpha,desca, Mdim, Ndim, W_dev)    ;  // imag contrib correction 
            commstr(Om_dev, &dP_dev[0]  , C1_dev, &neg_alpha,desca, Mdim, Ndim, W_dev)    ;  // imag contrib correction 
                                                                                             //        rmg::dcopy_driver(Nsq2, C_dev, ione, dP_dev ,ione)             ;  // dP=C     
            rmg::dcopy_driver(Nsq, C0_dev, ione, dP_dev ,ione)             ;  // dP=C     
            rmg::dcopy_driver(Nsq, C1_dev, ione, &dP_dev[Nsq] ,ione)             ;  // dP=C     
            rmg::daxpy_driver(Nsq2, rone,  dP_dev, ione, P1_dev, ione)       ;  // P1 =P1 +dP

            tstconv(dP_dev, &Nsq2, &thrs,&ierr,&err,&tConv, comm)  ;  // tstconv(dP,2*Nsq,N,thrs,ierr,err,tconv)
            if (iprint>0) rmg::printlog("ConvergTest: Niter  %d  errmax = %10.5e \n",  iter,err) ;
            if (abs(err) >  errmax)  errmax= abs(err)  ;

            iter ++ ;
        }

        rmg::dcopy_driver(Nsq2,  P1_dev, ione, P1, ione) ;   // dP =P0:  saves P0 into dP (for commutator )

        gpuFree(dP_dev);
        gpuFree(C1_dev);
        gpuFree(C0_dev);
        //gpuFree(C_dev);
        gpuFree(W_dev);
        gpuFree(P1_dev);
        gpuFree(P0_dev);
        gpuFree(Om_dev);

#endif
    }

    /*--- return  error max---- */
    *(p_errmax) = errmax ;
} 


//////////////////////////////////////////////////////////////////////////////////////
void commutp(std::complex<double> *P0, std::complex<double> *P1, std::complex<double> *Om, int *desca, int Mdim, int Ndim,
        double *p_thrs, int *p_maxiter, int *p_niter, double *p_errmax, int *p_iprint, MPI_Comm comm) 
{
    /*\
      !
      ! P = P0 +(dt/ih) [F,P]+ (dt/ih)^2/2! [F[F,P] + (dt/ih)^3/3! [F,F[F,P]] +... 
      !
      ! P(k)  <-- P(k-1) + 1/(ihk) *[Om,dP(k-1)] = 
      !         = P(k-1) +1/(hk)  *[Om,A ] -i/(hk) *[Om,S ]  
      !   
      !   where dP = S+i*A 
      !  
      !  Om = (Fav*dt)/I   from Magnus
      !
      \*/
    double   thrs    = *(p_thrs)       ;
    double   maxiter = *(p_maxiter)    ;
    int      iprint  = *(p_iprint)     ;

    int      Nsq = Mdim * Ndim         ;
    int      Nsq2 = 2*Mdim * Ndim      ;
    int      iter   = 1                ;
    bool     tConv  = false            ;
    double   errmax = 0.0e0            ;
    double   err    = 0.0e0            ;
    int      ierr   = 0                ;
    int      ione   = 1                ;

    // Needed for padding in C_dev

    std::complex<double>   rone(1.0, 0.0)           ;
    std::complex<double>   mone(-1.0,0.0)           ;

    /* ----  P1=P0 ,  dP=P0   ---- */
    std::complex<double> *C, *dP;
    MallocHostOrDevice((void **)&C,  Nsq * sizeof(std::complex<double>));
    MallocHostOrDevice((void **)&dP,  Nsq * sizeof(std::complex<double>));

    rmg::zcopy_driver(Nsq,  P0, ione, P1   ,ione) ;   // P1 =P0:  saves P0 into P1 (for updates) 
    rmg::zcopy_driver(Nsq,  P0, ione, dP   ,ione) ;   // dP =P0:  saves P0 into dP (for commutator )
    while ( iter  <= maxiter && tConv ==  false ) {

        std::complex<double> alpha(0.0, -1.0e0 /iter) ;
        std::complex<double> beta (0.0,0.0) ;
        if(ct.tddft_gpu && ct.tddft_tiledMM)
        {
            //C = Om * dP
            int Mglob = std::max(Mdim, Ndim);
            rmg::mgpu_zgemm_driver ("N", "N", Mglob, Mglob, Mglob, rone, Om, ione, ione, desca,
                    dP, ione, ione, desca, beta, C, ione, ione, desca);
            // dP = alpha * (C - C^dagger)
            rmg::zcopy_driver(Nsq,  C, ione, dP   ,ione) ;  
            rmg::zcommute_driver(alpha, Mglob, dP);

        }
        else
        {
            int Mglob = desca[3];
            if(ct.tddft_tiledMM) Mglob = std::max(Mdim, Ndim);
            rmg::mgpu_zgemm_driver ("N", "N", Mglob, Mglob, Mglob, alpha, Om, ione, ione, desca,
                    dP, ione, ione, desca, beta, C, ione, ione, desca);
            // C = -i * Om * dP

            rmg::mgpu_zgemm_driver ("N", "N", Mglob, Mglob, Mglob, alpha, dP, ione, ione, desca,
                    Om, ione, ione, desca, mone, C, ione, ione, desca);
            // C = -i (dP * OM - Om * dP)

            rmg::zcopy_driver(Nsq,  C, ione, dP   ,ione) ;  
        }
        //
        rmg::zaxpy_driver(Nsq, rone,  dP, ione, P1, ione)       ;  // P1 =P1 +dP
        tstconv((double *)dP, &Nsq2, &thrs,&ierr,&err,&tConv, comm)  ;  // tstconv(dP,2*Nsq,N,thrs,ierr,err,tconv)
        if (iprint>0) rmg::printlog("ConvergTest: Niter  %d  errmax = %10.5e \n",  iter,err) ;
        if (abs(err) >  errmax)  errmax= abs(err)  ;

        iter ++ ;
    }
    FreeHostOrDevice(C);
    FreeHostOrDevice(dP);
    /*--- return  error max---- */
    *(p_errmax) = errmax ;
} 
