  
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "blas.h"
#include <iostream>

using namespace std;

extern "C" { void diagev  (double*, double*, double*, int*, int*) ; }
extern "C" { void diagev_ (double*, double*, double*, int*, int*) ; }
extern "C" { void commutp_(double*, double*, double*, int*, double*, int*, int*, double*, int*) ;}
extern "C" { void fortfunc_(int *ii, double *ff); }
extern "C" { void eldyn_(int *N, double *S,double* F, double* Pn0,double* Pn1,int* Ieldyn,int* iprint); }



void diagev(double *, double *, double *, int *, int *);
void print_matrix(double *A, int N)  ;
void transpose(double *, double *, int);
void  commstr(double *A, double *B, double *C,double *p_alpha,int *p_Nbasis) ;
void commutp(double *P0, double *P1, double *Om, 
              int *p_Nbasis, double *p_thrs, int *p_maxiter, int *p_niter, double *p_errmax, int *p_iprint) ;

void eldyn_ort(int *p_N, double *F,double *Po0,double *Po1,int *p_Ieldyn,  double *thrs,int*maxiter,  double *errmax,int *niter , int *p_iprint) ;
void eldyn_nonort(int *p_N, double *S, double *F,double *Po0,double *Pn1,int *p_Ieldyn,  double *thrs,int*maxiter,  double *errmax,int *niter , int *p_iprint) ;

/////////////////////////////////////////////////////////////////////////
void print_matrix(double *A, int N) 
{ 
  int ij=0 ;
  for  (int i =0; i< N; i++ ) 
  {
     for  (int j =0; j< N; j++ ) { printf("  %14.8f",  A[ij++]); }
     printf("\n");
  }

}
/////////////////////////////////////////////////////////////////////////
void print_matrix2(double *A, int N) 
{ 
  
  int ij=0 ;
  printf(" **** real: \n");
  for  (int i =0; i< N; i++ ) 
  {
     for  (int j =0; j< N; j++ ) { printf("  %14.8f",  A[ij++]); }
     printf("\n");
  }
  printf("\n **** imag: \n");
  for  (int i =0; i< N; i++ ) 
  {
     for  (int j =0; j< N; j++ ) { printf("  %14.8f",  A[ij++]); }
     printf("\n");
  }


}
/////////////////////////////////////////////////////////////////////////
int  main()
{

 printf(" ECHO \n");

 int  Nbasis = 4;
 int Nbsq , Nbsq2 ;

 //double  fff =12.3456e1;
 Nbsq  =  Nbasis*Nbasis ;   
 Nbsq2 =  Nbsq*2 ;
   double *Hmatrix = new double[2*Nbsq ];
   double *Fn      = new double[2*Nbsq ];
//   double *Fn_dt   = new double[2*Nbsq ];
   double *Pn0     = new double[2*Nbsq ];
   double *Pn1     = new double[2*Nbsq ];
   double *S       = new double[  Nbsq ];
   double *S1      = new double[  Nbsq ];

   double zero = 0.0e0 ;
//   double  dt  = 1.0e0 ;


   for  (int i=0; i< Nbsq2  ; i++) { 
      Fn[i]  = zero ; 
      Pn0[i] = zero ; 
      Pn1[i] = zero ; 
   }
  printf(" === OK-3 \n");  fflush(stdout); 

/////////////////////////////////////////////////////////////////////////////////
/*
  Random hamiltonian:
 H =
  -0.5812069   0.8154454   0.8675571  -0.3967354
   0.8154454   0.0080746  -0.4762238  -0.1593811
   0.8675571  -0.4762238   0.3729276   0.3634973
  -0.3967354  -0.1593811   0.3634973   0.2870880

 Trandom density matrix: 
 P =

   0.057294  -0.194598   0.123567  -0.029571
  -0.194598   0.801742  -0.273484  -0.215150
   0.123567  -0.273484   0.418322  -0.391492
  -0.029571  -0.215150  -0.391492   0.722642
*/

Hmatrix[0]  = -0.5812069 ; Hmatrix[1]  =  0.8154454 ; Hmatrix[2]  =  0.8675571 ; Hmatrix[3]  = -0.3967354  ;
Hmatrix[4]  =  0.8154454 ; Hmatrix[5]  =  0.0080746 ; Hmatrix[6]  = -0.4762238 ; Hmatrix[7]  = -0.1593811  ;
Hmatrix[8]  =  0.8675571 ; Hmatrix[9]  = -0.4762238 ; Hmatrix[10] =  0.3729276 ; Hmatrix[11] =  0.3634973  ;
Hmatrix[12] = -0.3967354 ; Hmatrix[13] = -0.1593811 ; Hmatrix[14] =  0.3634973 ; Hmatrix[15] =  0.2870880  ;

/*  -- Pn0_re --- */
Pn0[0]  =   0.057294 ; Pn0[1]  =  -0.194598  ; Pn0[2]  =    0.123567  ; Pn0[3]  =   -0.029571  ; 
Pn0[4]  =  -0.194598 ; Pn0[5]  =   0.801742  ; Pn0[6]  =   -0.273484  ; Pn0[7]  =   -0.215150  ; 
Pn0[8]  =   0.123567 ; Pn0[9]  =  -0.273484  ; Pn0[10] =    0.418322  ; Pn0[11] =   -0.391492  ; 
Pn0[12] =  -0.029571 ; Pn0[13] =  -0.215150  ; Pn0[14] =   -0.391492  ; Pn0[15] =    0.722642  ; 

/*  -- Pn0_im --- */
Pn0[16] =    0.00000 ;  Pn0[17] =  0.11636  ; Pn0[18] = -0.24287  ; Pn0[19] =  0.06098  ; 
Pn0[20] =   -0.11636 ;  Pn0[21] =  0.00000  ; Pn0[22] =  0.12072  ; Pn0[23] = -0.42491  ; 
Pn0[24] =    0.24287 ;  Pn0[25] = -0.12072  ; Pn0[26] =  0.00000  ; Pn0[27] = -0.08322  ; 
Pn0[28] =   -0.06098 ;  Pn0[29] =  0.42491  ; Pn0[30] =  0.08322  ; Pn0[31] =  0.00000  ; 


S[0]  =   1.00000 ;  S[1]  = -0.28911 ;  S[2]  =  0.22935 ;  S[3]  =  0.14159 ;  
S[4]  =  -0.28911 ;  S[5]  =  1.00000 ;  S[6]  = -0.33655 ;  S[7]  = -0.20190 ;
S[8]  =   0.22935 ;  S[9]  = -0.33655 ;  S[10] =  1.00000 ;  S[11] =  0.14679 ;
S[12] =   0.14159 ;  S[13] = -0.20190 ;  S[14] =  0.14679 ;  S[15] =  1.00000 ;

S1[0]  =   1.00000 ;  S1[1]  = -0.00000 ;  S1[2]  =  0.00000 ;  S1[3]  =  0.00000 ;  
S1[4]  =   0.00000 ;  S1[5]  =  1.00000 ;  S1[6]  = -0.00000 ;  S1[7]  = -0.00000 ;
S1[8]  =   0.00000 ;  S1[9]  = -0.00000 ;  S1[10] =  1.00000 ;  S1[11] =  0.00000 ;
S1[12] =   0.00000 ;  S1[13] = -0.00000 ;  S1[14] =  0.00000 ;  S1[15] =  1.00000 ;



/////////////////////////////////////////////////////////////////////////////////
//
// Non-Orthogonal propagation:
//   

  printf("==================== INITIAL DATA ===========================\n");
  int ione = 1 ; 

  printf("  Hmatrix = \n");
  print_matrix(Hmatrix,Nbasis) ;

  //dt =10.0e0 ;
  dcopy_(&Nbsq, Hmatrix, &ione,  Fn,&ione) ;
  printf("  Fn*drt = \n");
  
  print_matrix(Fn,Nbasis) ;
 
  int iprint =0  ;

  printf("  S  = \n");
  print_matrix(S,Nbasis) ;


  printf("  Pn0  = \n");
  print_matrix2(Pn0,Nbasis);  


 // call commutp (Po0,Po1,Fo,N,thrs,maxiter,niter,errmax,iprint)
/****  call diagev (Po0,Po1,Fo,N,iprint)    ****/
//   diagev_(Po0,Po1, Fo, &Nbasis,&iprint) ;

  printf("==========================================================================\n");
  printf("\n\n  full FORTRAN   version of eldyn_non-orthogonal: \n\n\n"); 

  int Ieldyn = 1 ;  
  //  commutp: [in] :  thrs, maxiter   
  //           [out]:  errmax, niter
  //  diagp  : [in/out] :  <control=empty> 

  
 //  "eldyn.f90"  has hard-wired parameters :  thrs=1.0e-14 and maxiter=1000;
 //  maxiter and thrs  only affects  "eldyn_nonort.cpp"  and "eldyn_ort.cpp"
  int niter ;
  int     maxiter = 100 ;  
  double  thrs    = 1e-14  ;     
  double  errmax ;

  
 /* -------- Fortran  version: \n");  ---------*/ 

  Ieldyn =1 ;  iprint=4; 
  eldyn_(&Nbasis,S,Fn,Pn0,Pn1,&Ieldyn,&iprint)  ;
  printf("** Pn1 (after eldyn.f90)  :\n"); print_matrix2(Pn1,Nbasis);


  printf("==========================================================================\n");
  printf("\n\n  full C++ version of eldyn_non-orthogonal: \n\n\n"); 


  eldyn_nonort(&Nbasis, S,Fn,Pn0,Pn1,&Ieldyn, &thrs,&maxiter,  &errmax,&niter ,    &iprint) ;   
  printf("** Pn1 (after eldyn2.cpp)  :\n"); print_matrix2(Pn1,Nbasis);


//  printf(" Commutator test: \n"); 
//  double alpha=10.0e0 ;
//  commstr(Po0, Fo, &Po0[Nbsq], &alpha,&Nbasis)  ;
//

//  printf(" *** Po0 :\n");  print_matrix(Po0, Nbasis);
//  printf(" *** Fo  :\n");  print_matrix(Fo , Nbasis);
//  printf(" *** 10x[Po0,Fo] :\n");  print_matrix(&Po0[Nbsq], Nbasis);




   return 0;
   


}
