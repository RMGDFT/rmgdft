/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "rmgthreads.h"
#include "vhartree.h"
#include "packfuncs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "Functional.h"
#include "Solvers.h"
#include "RmgParallelFft.h"

#include "blas_driver.h"
#include "prototypes_tddft.h"
#include "RmgException.h"


///////////////////////////////////////
void  magnus( double *H0, double *H1, double p_time_step , double *Hdt, int ldim){
    /*  
        This is a first order Magnus  expansion, wchih is equivalent to   mid-point propagator
        F*dt =   Integrate_0^t   H(t) dt   
        ~= 1/2*(H0 +H1) *dt

        It works for real  Hamiltonian matrices only.  Complex are required for bybrid functions with exact exchange

     */  

    double  dt     =    p_time_step  ;   // time step
    int     Nbasis =    ldim         ;     // leaading dimension of matrix

    int     Nsq =    Nbasis*Nbasis    ; 
    //double  dminus_one  = -1.0e0      ; 
    double  one         =  1.0e0      ; 
    double  half_dt     =  0.5e0*dt   ;
    int    ione = 1 ;
    //printf("magnus,  dt, ldim =  %f   %d  \n", dt,Nbasis ) ;

    dcopy_driver ( Nsq ,   H0      ,     ione , Hdt ,  ione) ;
    daxpy_driver ( Nsq ,  one      , H1, ione , Hdt ,  ione) ;   
    dscal_driver ( Nsq ,  half_dt  ,            Hdt ,  ione) ;


}
///////////////////////////////////////
void tst_conv_matrix  (double * p_err , int * p_ij_err ,   double *H0, double *H1,  int ldim) {  
    /* 
       Test convergemce: calculate  infty error for convergence:   || H1 - H0||_infty
       Here error is  L_infty norm of a difference matrix dH =H1-H0:
     *p_err    :  [out]  returns absolute value of the largest matrix element of  dH  where dH = H1-H0 
     *p_ij_err :  [out]  returns position  of err in the  matrix dH
     H0, H1    :  [in]   tested matrices
ldim      :  [in]   leading dimension of  H0, H1  (= Nbasis)
     */


    int    Nbasis     =   ldim         ;  
    int    Nsq        =   Nbasis*Nbasis ;  

    double err        =  fabs( H1[0] - H0[0]) ; 
    int    ij_err     =  0                    ;  // position/ location  in matrix

    for (int i =1 ;  i < Nsq ; i++ )  {
        double  tst = fabs(H1[i]-H0[i])  ;
        if (  tst  > err)   { err = tst ; ij_err = i ; } 
    } 

    //  pass back the max  error and its location:
    * p_err    =  err   ;
    * p_ij_err = ij_err ; 

} 
///////////////////////////////////////
void extrapolate_Hmatrix   (double  *Hm1, double *H0, double *H1,  int ldim) {  
    /* 
       Linear extrapolation of  H(1) from H(0) and  H(-1):
       H(1) =   H(0) + dH   =  H(0) + ( H(0)-H(-1) )  =  2*H(0) - H(-1) 
       H1   =  2*H0 - Hm1

       Hm1,H0  : [in]
H1      : [out]
     */

    int    Nbasis     =   ldim         ;  
    int    Nsq        =   Nbasis*Nbasis ;  

    int    ione       = 1 ;
    double  neg_one   = -1.0e0   ;
    double  two       =  2.0e0   ;

    dcopy_driver ( Nsq ,           Hm1 , ione , H1 ,  ione) ;    //  
    dscal_driver ( Nsq , neg_one , H1  , ione )              ;    //  H1 = -H(-1)  
    daxpy_driver ( Nsq , two     , H0  , ione , H1  , ione) ;    //  H1 =  2*H0 + H1 =  2*H0 -Hm1 

}
