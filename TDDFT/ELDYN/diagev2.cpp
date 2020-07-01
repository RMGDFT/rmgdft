 //  diagev2_(Po0,Po1, Fo, &Nbasis,&iprint) ;
  
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


////////////////////////////////////////////////////////////
void set_vector( double scalar, double  *A, int N) 
{
   //int ij=0 ;
  for (int i=0; i<N ; i++) { A[i] =scalar ; }
}


////////////////////////////////////////////////////////////

void diagev( double *Po0, double *Po1, double *Fdt, int * p_Nbasis,  int * p_iprint )
{
  int Nbasis = *(p_Nbasis) ;
  //int iprint = *(p_iprint) ;

  int Nbsq  =   Nbasis*Nbasis  ;
  int Nbsq2 = 2*Nbasis*Nbasis  ;

  int ione  = 1 ;
  double one =1.0e0 , minus_one=-1.0e0 ,  zero = 0.0e0 ;

  //allocate( C(N,N), eig(N), work(5*N*N) )
  double *eig   = new double[ Nbasis ] ;
  double *C     = new double[ Nbsq2  ] ; 
  double *work  = new double[ 5*Nbsq ] ;
  int     lwork =             5*Nbsq   ;
  int     info  ;

  dcopy_(&Nbsq2, Po0, &ione,Po1, &ione) ;
  dcopy_(&Nbsq2, Fdt, &ione,C  , &ione) ;

  set_vector(zero,Po1,Nbsq2) ; //print_matrix2(Po1,Nbasis);

    char jobz='V' ; 
    char uplo='U' ;
    
    //dsyev_ ( 'Vectors', 'Upper', &Nbasis, C, &Nbasis, eig, work, &lwork, &info ) ;
    dsyev_ ( &jobz,  &uplo , &Nbasis, C, &Nbasis, eig, work, &lwork, &info ) ;

/*\
!
!  time evolution operator
! U  = CC * exp(   F*dt/i*hbar) * CC' 
!    = CC * exp( -i*Fdt/hbar  ) * CC'  = C*expRe*C'  +  i* C*expIm*C'
!
\*/  

  double *Ure  = new double[ Nbsq  ] ;    
  double *Uim  = new double[ Nbsq  ] ;    
  double *Wre  = new double[ Nbsq  ] ;    
  double *Wim  = new double[ Nbsq  ] ;    

  set_vector(zero, Ure, Nbsq) ;
  set_vector(zero, Uim, Nbsq) ;
  set_vector(zero, Wre, Nbsq) ;
  set_vector(zero, Wim, Nbsq) ;

/*\
!
!  W=C*exp  
!
\*/

 for (int i=0;i<Nbasis ; i++)  
 {
    double expRe = cos(-eig[i]) ;
    double expIm = sin(-eig[i]) ;
    int j = i*Nbasis  ;

    daxpy_(&Nbasis,&expRe,&C[j], &ione, &Wre[j], &ione) ;
    daxpy_(&Nbasis,&expIm,&C[j], &ione, &Wim[j], &ione) ;
 }  
/*\
!
!  U =W*C^T
!
\*/
  int N = Nbasis ; 
   char trN='N' ;  // transpose ='not'
   char trT='T' ;  //  transpose='True'
   
   dgemm_(&trN, &trT, &N,&N,&N, &one, Wre, &N, C, &N, &zero, Ure,&N) ;  
   dgemm_(&trN, &trT, &N,&N,&N, &one, Wim, &N, C, &N, &zero, Uim,&N) ;  
 
/*\
!
! now propagate:  P1 = U* P0 *U' 
!       where U=(Ur+i*Ui)  and  U'=(Ur   -i*Ui^T)
!
!   P1 = (Ur+i*Ui)*(Pr +i*Pi) *( Ur -i*Ui^T) 
!
! dgem: C := alpha*op( A )*op( B ) + beta*C, 
!  DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 
!
!   double * p_P0r = &Po0[0]    ;
!   double * p_P0i = &Po0[Nbsq] ;
!   double * p_P1r = &Po1[0]    ;
!   double * p_P1i = &Po1[Nbsq] ;
\*/

  dgemm_(&trN,&trN,&N,&N,&N,       &one,Ure,&N,&Po0[0]   ,&N,&zero, Wre,      &N);    //  Wre  =   Ur*P0r 
  dgemm_(&trN,&trN,&N,&N,&N,       &one,Wre,&N, Ure      ,&N,&zero, Po1,      &N);    //  P1r  =   UrP0r*Ur 
  dgemm_(&trN,&trT,&N,&N,&N, &minus_one,Wre,&N, Uim      ,&N,&zero,&Po1[Nbsq],&N);    //  P1i  = -i*UrP0r*Ur 
                                                                    
  dgemm_(&trN,&trN,&N,&N,&N,       &one,Uim,&N,&Po0[0]   ,&N,&zero, Wim,      &N);    //  Wim  =  i*Ui*Pr 
  dgemm_(&trN,&trN,&N,&N,&N,       &one,Wim,&N, Ure      ,&N,&one, &Po1[Nbsq],&N);    //  P1i +=  i*UiPr*Ur
  dgemm_(&trN,&trT,&N,&N,&N,       &one,Wim,&N, Uim      ,&N,&one,  Po1,      &N);    //  P1r +=  i*UiPr*(-i*Ui^t) = +UiPr*UiT
                                                                    
  dgemm_(&trN,&trN,&N,&N,&N,       &one,Ure,&N,&Po0[Nbsq],&N,&zero, Wim,      &N);    //  Wim  =  i*Ur*Pi   
  dgemm_(&trN,&trN,&N,&N,&N,       &one,Wim,&N, Ure      ,&N,&one ,&Po1[Nbsq],&N);    //  P1i +=  i*UrPi*Ur  
  dgemm_(&trN,&trT,&N,&N,&N,       &one,Wim,&N, Uim      ,&N,&one , Po1,      &N);    //  P1r +=  i*UrPi*(-i*Ui^t) = +UrPi*UiT
                                                                    
  dgemm_(&trN,&trN,&N,&N,&N, &minus_one,Uim,&N,&Po0[Nbsq],&N,&zero, Wre,      &N);    //  Wre  = (i*i)*Ui*Pi       =   -UiPi
  dgemm_(&trN,&trN,&N,&N,&N,       &one,Wre,&N, Ure      ,&N,&one , Po1,      &N);    //  P1r += (-UiPi)*Ur 
  dgemm_(&trN,&trT,&N,&N,&N, &minus_one,Wre,&N, Uim      ,&N,&one ,&Po1[Nbsq],&N);    //  P1i += (-UiPi)*(-i*Ui^t)
  
 // printf("*** Final Po1 :\n"); print_matrix2(Po1,Nbasis);

}


