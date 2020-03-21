#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "blas.h"
#include <iostream>




void commutp(double *P0, double *P1, double *Om,
                int *p_Nbasis, double *p_thrs, int *p_maxiter, int *p_niter, double *p_errmax, int *p_iprint) ;
  
void diagev(double *, double *, double *, int *, int *);
void dochol(double *S,  double *Si, double *U, double *Ui,  int*N) ;  
void ontrans(double *U,double *Ui, double*Ain, double*Aout,int*N,int*ifon) ;
void fill_LT(double *A, int *p_N, int *p_lt) ;
void print_matrix(double*, int);
void print_matrix2(double*, int);
void set_vector(double zero, double *Wre, int Nbsq) ;
void clean_LT(double *A, int *p_N, int *p_lt) ;

void ontrans(double*U, double*Ui, double*Ain, double *Aout, int*p_N, int*p_ifon) ;

//void eldyn_ort(int *p_N, double *F,double *Po0,double *Po1,int *p_Ieldyn,  double *thrs,int*maxiter,  double *errmax,int *niter , int *p_iprint) ;
    
////////////////////////////////////////////////

void eldyn_nonort(int *p_Nbasis, double *S,  double *Omega,double *Pn0,double *Pn1,int *p_Ieldyn, double *p_thrs, int *p_maxiter,  double *p_errmax, int *p_niter, int*p_iprint) 
{
/*
  Omega  [in]  : 
  Po0    [in]  : assumed in orhotgonal basis set
  Po1    [out] : assumed in orthogonal basis set
  Ieldyn [in]  : 1= use BCH expansion, 2=  use diagonalization , 0= use both  BCH and diag 
*/
 int     Nbasis   = *(p_Nbasis)       ;   //  Nbasis    
 int     Ieldyn   = *(p_Ieldyn)  ;
 int     maxiter  = *(p_maxiter) ;
 double  thrs     = *(p_thrs)    ;
 int     iprint   = *(p_iprint)  ;

 bool  tCommutp                  ;
 bool  tDiagp                    ;

 int    niter  ;
 double errmax ;

 int N = Nbasis ;
 int Nsq  = N*N ;
 int Nsq2 = 2*Nsq ;



 /*-- select propagator --*/ 
   if (Ieldyn == 0 )       {
       tCommutp = true  ;
       tDiagp   = true  ;
      //printf("Ieldyn = 0 \n");

   } else if (Ieldyn == 1) {
       tCommutp = true  ;
       tDiagp   = false ;
      //printf("Ieldyn = 1 \n");

   } else if (Ieldyn == 2) {
       tCommutp = false ;
       tDiagp   = true  ;
      //printf("Ieldyn = 2 \n");

   } else  {
      printf("Incorrect ieldyn =%d. Should  be 0, 1 or 2. Exiting.\n", Ieldyn);
      return ;
   }


/*----- get space for orthogonalizatio ----*/
  double *Si = new  double[Nsq ] ;
  double *U  = new  double[Nsq ] ;
  double *Ui = new  double[Nsq ] ;
  double *Omega_o   = new  double[Nsq2 ] ;
  double *Po0       = new  double[Nsq2 ] ;
  double *Po1       = new  double[Nsq2 ] ;

// set o/n transformation flags:
   int  iPo2Pn = 1  ; 
   int  iPn2Po = 2  ; 
   int  iFo2Fn = 3  ; 
   int  iFn2Fo = 4  ;

/*----- get transformation matrices ----*/
    dochol(S,Si,U,Ui,&N)  ;  
    ontrans(U,Ui, Omega   ,  Omega_o ,&N,&iFn2Fo) ;
    ontrans(U,Ui,&Pn0[0]  , &Po0[0]  ,&N,&iPn2Po) ;
    ontrans(U,Ui,&Pn0[Nsq], &Po0[Nsq],&N,&iPn2Po) ;
/* 
    printf("*** U  (c++) : \n")       ; print_matrix(U ,     N) ;
    printf("*** Ui (c++) : \n")       ; print_matrix(Ui,     N) ;
    printf("*** Si (c++) : \n")       ; print_matrix(Si,     N) ;
    printf("*** F/Omega_o (c++):x \n"); print_matrix(Omega_o,N) ;
    printf("*** Po0 (c++) :\n")       ; print_matrix2(Po0,   N) ;
*/

 /*-----  run it in orthogonal  basis set -----*/
  if (tCommutp) { 
      commutp (Po0,Po1,Omega_o,&Nbasis,&thrs,&maxiter,&niter,&errmax,&iprint) ;
  } 
  if (tDiagp)   {
      diagev(Po0, Po1,Omega_o ,&Nbasis,&iprint);
  } 

/*  --  now back n/o transformation: Po->Pn  --- */
    ontrans(U,Ui,&Po1[0]  , &Pn1[0]  , &N,&iPo2Pn) ;
    ontrans(U,Ui,&Po1[Nsq], &Pn1[Nsq], &N,&iPo2Pn) ;


    //printf("*** Po1 (c++) :\n")       ; print_matrix2(Po1,   N) ;
    //printf("*** Pn1 (c++) :\n")       ; print_matrix2(Pn1,   N) ;
    //printf(" ----- leaving eldyn_nonort --- \n");
}
