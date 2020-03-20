  
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "blas.h"
#include <iostream>
//#include "lapacke.h"

void print_matrix(double *, int ) ;
void print_matrix2(double *, int ) ;

/////////////////////////////////////////////////////////////////////////
void transpose( double *A,  double *B, int Nbasis) 
{
  //printf(" transpose  %d \n", Nbasis) ;
  int ij=0  ;
  for (int ix=0; ix< Nbasis;ix++) {
     for (int iy=0; iy<Nbasis;  iy++) {
         B[iy*Nbasis +ix] =A[ij++]  ;
     }   
  }
  //printf(" *** A :     \n"); print_matrix(A,Nbasis);
  //printf(" *** B= A^T :\n"); print_matrix(B,Nbasis);
}

/////////////////////////////////////////////////////////////////////////

void  commstr(double *A, double *B, double *C,double *p_alpha,int *p_Nbasis)
{
  /*\
  ! commutator  of two symmetric matrices:
  !   [S1,S2] =S1*S2-*S2*S1= S1*S2 -transpose(S1*S2)
  \*/
 double  alpha  = *(p_alpha)  ;
 int     Nb     = *(p_Nbasis) ;
 double  beta   = 0.0e0       ;
 int     Nbsq   = Nb *Nb      ;
  
 double *W   = new double[ Nbsq ] ;
 char trN='N' ;  // transpose ='not
 // call dgemm('N','N',Nb,Nb,Nb,alpha,A,Nb,B,Nb,beta,C,Nb)
 dgemm_(&trN,&trN,&Nb,&Nb,&Nb,&alpha,A,&Nb,B,&Nb,&beta,C,&Nb) ;
 transpose(C,W,Nb) ;

//SUBROUTINE DAXPY    ( n, alpha, x, incx, y, incy ) :: //  y  <--  alpha*x + y
 double minus_one = -1.0e0 ;
 int    ione      =  1     ;
 
 daxpy_(&Nbsq,  &minus_one , W, &ione, C, &ione) ;

}

/////////////////////////////////////////////////////////////////////////
void  commatr(double *A, double *B, double *C,double *p_alpha,int *p_Nbasis)
{
  /*\
  ! commutator of two anti-symmetric matrices:
  !  [S,A] = S*A-A*S =  S*A +transpose(S*A)
  \*/
  double  alpha  = *(p_alpha)  ;
  int     Nb     = *(p_Nbasis) ;
  double  beta   = 0.0e0       ;
  int     Nbsq   = Nb *Nb      ;
   
  double *W   = new double[ Nbsq ] ;
  char trN='N' ;  // transpose ='not
  // call dgemm('N','N',Nb,Nb,Nb,alpha,A,Nb,B,Nb,beta,C,Nb)
  dgemm_(&trN,&trN,&Nb,&Nb,&Nb,&alpha,A,&Nb,B,&Nb,&beta,C,&Nb) ;
  transpose(C,W,Nb) ;

  double one = 1.0e0 ;
  int    ione      =  1     ;
  daxpy_(&Nbsq,  &one , W, &ione, C, &ione) ;
}

//////////////////////////////////////////////////////////////////////////////////////

void  tstconv(double *C,int *p_M,int *p_N, double *p_thrs,int *p_ierr, double *p_err, bool *p_tconv) 
{

   int     M     = *(p_M)    ;  //  [in]  :  total  size of matrix (2*Nbasis*Nbasis)
   int     N     = *(p_N)    ;  //  [in]  :  Nbasis (not used)
   double  thrs  = *(p_thrs) ;  //  [in]  :  convergence threshold
   
   double  err   = abs(C[0]) ;  //  [out] :   error= abs of max element in the matrix
   int    ierr   =   0       ;  //  [out] :   location of err in matrix/vector
   bool   tconv  =  false    ;  //  [out] :   if converged ?  true or false?

   //err = abs(C[0]); 
   for (int i=0; i <M ;i++) {
     double err_tmp = abs(C[i]) ; 
     if (err_tmp > err) {
         err  = err_tmp ;
         ierr = i       ;
     }
   }

   if (err < thrs)  tconv = true ;
   /*-- return values **/
   *(p_err )  =  err  ;
   *(p_ierr)  =  ierr ;
   *(p_tconv) =  tconv ;
}
//////////////////////////////////////////////////////////////////////////////////////
void commutp(double *P0, double *P1, double *Om,
              int *p_Nbasis, double *p_thrs, int *p_maxiter, int *p_niter, double *p_errmax, int *p_iprint) 
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
  int      Nbasis  = *(p_Nbasis)        ;
  double   thrs    = *(p_thrs)          ;
  double   maxiter = *(p_maxiter)       ;
  int      niter   = *(p_niter)         ;
  int      iprint  = *(p_iprint)        ;
                                       
  int            N = Nbasis             ;
  int          Nsq =   N*N              ;
  int         Nsq2 = 2*Nsq              ;

  double      *C   = new double[ Nsq2 ] ;
  double      *dP  = new double[ Nsq2 ] ;

  int       iter   = 1                  ;
  bool      tConv  = false              ;
  double    errmax = 0.0e0              ;
  double    err    = 0.0e0              ;
  int       ierr   = 0                  ;
  int       ione   = 1                  ;
  double    one    = 1.0e0              ;

/* ----  P1=P0 ,  dP=P0   ---- */
 dcopy_(&Nsq2,  P0, &ione, P1   ,&ione) ;   // P1 =P0:  saves P0 into P1 (for updates) 
 dcopy_(&Nsq2,  P0, &ione, dP   ,&ione) ;   // dP =P0:  saves P0 into dP (for commutator )

 
 while ( iter  <= maxiter && tConv ==  false ) {
    double alpha     =  1.0e0 /iter ;
    double neg_alpha = -alpha       ;

    commatr(Om, &dP[Nsq], &C[0]  , &alpha,    &N)    ;  // real contrib correction  
    commstr(Om, &dP[0]  , &C[Nsq], &neg_alpha,&N)    ;  // imag contrib correction 

    dcopy_(&Nsq2, C, &ione, dP   ,&ione)             ;  // dP=C     
    daxpy_(&Nsq2, &one,  dP, &ione, P1, &ione)       ;  // P1 =P1 +dP
    tstconv(dP, &Nsq2 ,&N, &thrs,&ierr,&err,&tConv)  ;  // tstconv(dP,2*Nsq,N,thrs,ierr,err,tconv)
    if (iprint>0) printf("ConvergTest: Niter  %d  errmax = %10.5e \n",  iter,err) ;
    if (abs(err) >  errmax)  errmax= abs(err)  ;

    iter ++ ;
 }

/*--- return  error max---- */
  *(p_errmax) = errmax ;
} 
