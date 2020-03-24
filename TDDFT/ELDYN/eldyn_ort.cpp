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
    
////////////////////////////////////////////////

void eldyn_ort(int *p_Nbasis, double *Omega,double *Po0,double *Po1,int *p_Ieldyn, double *p_thrs, int *p_maxiter,  double *p_errmax, int *p_niter, int*p_iprint) 
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

 /*-- select propagator --*/ 
   if (Ieldyn == 0 )       {
       tCommutp = true  ;
       tDiagp   = true  ;
      // printf("Ieldyn = 0 \n");

   } else if (Ieldyn == 1) {
       tCommutp = true  ;
       tDiagp   = false ;
     // printf("Ieldyn = 1 \n");

   } else if (Ieldyn == 2) {
       tCommutp = false ;
       tDiagp   = true  ;
     // printf("Ieldyn = 2 \n");

   } else  {
      printf("Incorrect ieldyn =%d. Should  be 0, 1 or 2. Exiting.\n", Ieldyn);
      return ;
   }

 
 /*----- now run it  -----*/
  if (tCommutp) { 
     //commutp (Po0,Po1,Fo,N,thrs,maxiter,niter,errmax,iprint)
      commutp (Po0,Po1,Omega,&Nbasis,&thrs,&maxiter,&niter,&errmax,&iprint) ;
  } 
  if (tDiagp)   {
     // diagev2(Po0, Po1,Fo,&Nbasis,&iprint); 
    diagev(Po0, Po1,Omega ,&Nbasis,&iprint);

  } 



}
