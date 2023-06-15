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
    int sizeofmatrix, item, item1, item2, lwork;
    int nproc, myrow, mycol, icrow, iccol;
    int izero = 0, ione = 1, itwo = 2, nb, nn, NN;
    int nprow = pct.scalapack_nprow, npcol = pct.scalapack_npcol, npes = pct.grid_npes;
    int locr, qrmem, sizemqrleft, ldc, mpc0, nqc0, nrc;
    int NB;

    nb = ct.scalapack_block_factor;
    NB = ct.scalapack_block_factor;
    nn = ct.num_states;
    NN = ct.num_states;
    
    int S0_BASIS = (get_PX0_GRID()+2) * (get_PY0_GRID()+2) * (get_PZ0_GRID()+2);

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
    my_malloc_init( vtot_global, get_NX_GRID() * get_NY_GRID() * get_NZ_GRID(), double );
    my_malloc_init( wave_global, get_NX_GRID() * get_NY_GRID() * get_NZ_GRID(), double );
    rho_global = vtot_global;
    my_malloc_init( rho_old, get_FP0_BASIS() , double );

    my_malloc_init( sg_res, S0_BASIS, double );

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
    /*added by shuchun wang, it is used to calculate partial_omega/partial_R 
       for nonlocal force */
    my_malloc_init( mat_Omega, sizeofmatrix, double );

    /*  allocate memory for other uses  */

    item = ct.num_kpts * ct.num_states;
    item1 = item;
    item2 = 2 * get_P0_BASIS() + S0_BASIS + item1;
    item = rmg_max(13 * S0_BASIS, item2);

    nproc = pct.scalapack_nprow * pct.scalapack_npcol;
    locr = ((ct.num_states / NB + 1) / nproc + 1) * NB + NB;
    lwork = 10 * (locr * 5 + NB);
    item1 = rmg_max(lwork, item);

    nprow = pct.scalapack_nprow;
    npcol = pct.scalapack_npcol;
    myrow = pct.scalapack_myrow;
    mycol = pct.scalapack_mycol;

    icrow = indxg2p(&itwo, &nb, &myrow, &izero, &nprow);
    iccol = indxg2p(&ione, &nb, &mycol, &izero, &npcol);
    mpc0 = numroc(&nn, &nb, &myrow, &icrow, &nprow);
    nqc0 = numroc(&nn, &nb, &mycol, &iccol, &npcol);
    nrc = numroc(&nn, &nb, &myrow, &izero, &npes);
    ldc = rmg_max(1, nrc);
    sizemqrleft = rmg_max((NB * (NB - 1)) / 2, (nqc0 + mpc0) * NB) + NB * NB;
    sizemqrleft *= 2;
    qrmem = 2 * NN - 2;
    lwork = 5 * NN + NN * ldc + rmg_max(sizemqrleft, qrmem) + 1;

    item = rmg_max(lwork, item1);
    my_malloc_init( work_memory, item, double );



}                               /* end allocate_matrix */
