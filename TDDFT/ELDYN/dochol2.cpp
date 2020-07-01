  
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

///////////////////////////////////////////////////////////////
void clean_LT(double *A, int *p_N, int *p_lt)
{
   int N  = *(p_N) ;
   int lt = *(p_lt);

   double zero =0.0e0 ;

   if (lt == 1 ) {
      for (int i=0; i< N; i++) {
          for (int j=0;j<i;j++){
             int ji=j*N+i ; 
              A[ji] =zero ;
          }
      }
   } else if (lt ==2) {
      for (int i=0; i< N; i++) {
          for (int j=0;j<i;j++){
             int ij=i*N+j ; 
              A[ij] =zero ;
          }
      }
   }
}

///////////////////////////////////////////////////////////////
void fill_LT(double *A, int *p_N, int *p_lt) 
{
/*
subroutine fillLT(A,N,lt)
!  lt=1 : fill L part by copying from U
!  lt=2 : fill U part by copying from L
*/
  int N  = *(p_N) ;
  int lt = *(p_lt);

  int upper_to_lower = 1 ;  // copy upper to lower diagonal 
  int lower_to_upper = 2 ;  // copy lower to upper diagonal 

  //double zero  =0.0e0 ;

  if (lt == upper_to_lower ) {
     for (int i=0;i< N;i++) {
        for (int j=0;j<i; j++) {
           int ij = i*N+j;
           int ji = j*N+i;
           A[ji] =A [ij];
        }
     } 
  } else if (lt == lower_to_upper) {
     for (int i=0;i< N;i++) {
        for (int j=0;j<i; j++) {
           int ij = i*N+j;
           int ji = j*N+i;
           A[ij] =A [ji];
        }
     } 
  }
}


///////////////////////////////////////////////////////////////
 void  ontrans(double*U, double*Ui, double*Ain, double *Aout, int*p_N, int*p_ifon)
{
/*\  subroutine  ontrans (U,Ui, Ain, Aout,N,ifon)
 !-------------------------------------------------
 !    Sn = L*U        ;  Sni = Ui*Li
 !    So = Li*Sn*Ui   ; 
 !    Sn = L* I* U    ; 
 !    Fn = L*Fo* U    ;  Fo = Li*Fn*Ui  
 !    Pn = Ui*Po*Li   ;  Po = U *Pn*L
 !------------------------------------------------
 ! ifon=1 :   Po->  Pn       :     Ui * Po * (Ui^t)
 !      2 :   Pn->  Po       :     U  * Pn * (U^t)
 !      3 :   Fo->  Fn       :  (U^t) * Fo *  U
 !      4 :   Fn->  Fo       : (Ui^t) * Fn *  Ui
 !--------------------------------------------------
\*/ 

  int N       = *(p_N)  ;
  int ifon    = *(p_ifon)    ;

  double zero = 0.0e0 ;
  double one  = 1.0e0 ;


  int   Nsq=N*N  ;
  double * W  = new double [Nsq] ;

  char  fN = 'N' ;
  char  fT = 'T' ;

  if (ifon == 1) {                //    !  Po->Pn  :  Pn= Ui * Po * (Ui^t)
          dgemm_(&fN,&fN,&N,&N,&N,&one,Ui,&N,Ain,&N,&zero,W   ,&N)  ;
          dgemm_(&fN,&fT,&N,&N,&N,&one,W ,&N,Ui ,&N,&zero,Aout,&N)  ;
  } else if (ifon == 2) {        //    !  Pn->Po  :  Po= U  * Pn * (U^t)
          dgemm_(&fN,&fN,&N,&N,&N,&one,U ,&N,Ain,&N,&zero,W   ,&N)  ;
          dgemm_(&fN,&fT,&N,&N,&N,&one,W ,&N,U  ,&N,&zero,Aout,&N)  ;

  } else if (ifon == 3) {          //    !  Fo->Fn  :  Fn= (U^t) * Fo *  U 
          dgemm_(&fT,&fN,&N,&N,&N,&one,U ,&N,Ain,&N,&zero,W   ,&N)  ;
          dgemm_(&fN,&fN,&N,&N,&N,&one,W ,&N,U  ,&N,&zero,Aout,&N)  ;

  } else if (ifon == 4)  {         //    !  Fn->Fo  :  Fo= (Ui^t) * Fn *  Ui
          dgemm_(&fT,&fN,&N,&N,&N,&one,Ui,&N,Ain,&N,&zero,W   ,&N)  ;
          dgemm_(&fN,&fN,&N,&N,&N,&one,W ,&N,Ui ,&N,&zero,Aout,&N)  ;
  } else  {
     printf(" Wrong  o/n transformation\n") ;
     exit( 1) ; 
  }

}

///////////////////////////////////////////////////////////////
 void  dochol(double *S,double *Si,double *U,double *Ui, int *p_N)
{
  int N = *(p_N); 
  int Nsq = N*N ; 
  int lt = 1 ; 
  int info  ; 

  //double *W =  new double[Nsq] ; 
  int ione =1 ;
  char  fU ='U';
  char  fN ='N';

/*---    do cholesky and  find inverse of S  --- */
  dcopy_(&Nsq,S,&ione,Si,&ione) ;
  dpotrf_(&fU,&N,Si,&N,&info)     ;    //!  cholesky 
  if (info != 0) { printf(" 'Chol-1-err ;stop;\n"); }
  dpotri_(&fU ,&N,Si,&N,&info)    ;    //  inverse  of S using  cholesky 
  if (info != 0)  { printf(" 'Chol-2-err';stop; \n"); }
  fill_LT(Si,&N,&lt);

//
//--- now get  transf.  matrices
//
  
  dcopy_(&Nsq,S,&ione,U,&ione) ;
  dpotrf_(&fU,&N,U,&N,&info)    ;  // !  cholesky 
  if (info != 0) { printf(" 'Chol-1-err ;stop;\n"); }
  dcopy_(&Nsq,U,&ione,Ui,&ione) ;
  dtrtri_(&fU,&fN,&N,Ui,&N,&info) ;
  if (info != 0) { printf(" 'Chol-4-err';stop \n"); }
  clean_LT(U,&N,&lt) ;
  clean_LT(Ui,&N,&lt);
}
