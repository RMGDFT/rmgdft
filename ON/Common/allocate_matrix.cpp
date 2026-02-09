/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/allocate_matrix.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
*                       Jerzy Bernholc
 * FUNCTION
 *   void allocate_matrix()   
 *   allocate memory for matrixs
 * INPUTS
 *   nothing
 * OUTPUT
 *   nothing
 * PARENTS
 *   run.c
 * CHILDREN
 
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"


void allocate_matrix()
{
    int sizeofmatrix;

    my_malloc_init( rho, get_FP0_BASIS() *2, double );
    rho_oppo = rho + get_FP0_BASIS() ;
    my_malloc_init( rho_tot, get_FP0_BASIS(), double );
    my_malloc_init( rhoc, get_FP0_BASIS(), double );
    my_malloc_init( vh, 3*get_FP0_BASIS(), double );
    vxc = vh + get_FP0_BASIS();
    my_malloc_init( vnuc, get_FP0_BASIS(), double );
    my_malloc_init( vtot, get_FP0_BASIS(), double );
    my_malloc_init( vtot_c, get_P0_BASIS(), double ); /*shuchun add */
    my_malloc_init( rhocore, get_FP0_BASIS(), double );
    my_malloc_init( eig_rho, get_FP0_BASIS(), double );
    my_malloc_init( rho_old, get_FP0_BASIS() , double );

    sizeofmatrix = MXLLDA * MXLCOL;

    if(!ct.is_gamma)
        sizeofmatrix *= 2.0;


    my_malloc_init( statearray, sizeofmatrix, double );
    my_malloc_init( l_s, sizeofmatrix, double );
    my_malloc_init( matB, sizeofmatrix, double );
    my_malloc_init( mat_X, sizeofmatrix, double );
    my_malloc_init( Hij, sizeofmatrix, double );
    my_malloc_init( work_dis, sizeofmatrix, double );
    my_malloc_init( work_dis2, sizeofmatrix, double );
    my_malloc_init( zz_dis, sizeofmatrix, double );
    my_malloc_init( cc_dis, sizeofmatrix, double ); //transpose of zz_dis
    my_malloc_init( gamma_dis, sizeofmatrix, double );
    my_malloc_init( uu_dis, sizeofmatrix, double );


}                               /* end allocate_matrix */
