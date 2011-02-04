/************************** SVN Revision Information **************************
 **    $Id: allocate_matrix_soft.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/allocate_matrix_soft.c *****
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
 *   
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "md.h"


void allocate_matrix_soft ()
{
    int ispin, sizeofmatrix, item, item1, item2, lwork;
    int ictxt, nproc, myrow, mycol, icrow, iccol;
    int izero = 0, ione = 1, itwo = 2, nb = NB, nn = NN;
    int nprow = pct.nprow, npcol = pct.npcol, npes = NPES;
    int locr, qrmem, sizemqrleft, ldc, mpc0, nqc0, nrc;

    ispin = ct.spin + 1;

    my_malloc_init( rho, FP0_BASIS * ispin, REAL );
    my_malloc_init( rhoc, FP0_BASIS, REAL );
    my_malloc_init( vh, FP0_BASIS, REAL );
    my_malloc_init( vnuc, FP0_BASIS, REAL );
    my_malloc_init( vxc, FP0_BASIS * ispin, REAL );
    my_malloc_init( vtot, FP0_BASIS, REAL );
    my_malloc_init( vtot_c, P0_BASIS, REAL );
    my_malloc_init( rhocore, FP0_BASIS, REAL );
/*  my_malloc_init( nlarray1, P0_BASIS/MAX_FUNC_PE, REAL );*/
/*  my_malloc_init( vtot_global, NX_GRID * NY_GRID * NZ_GRID, REAL );
 * rho_global = vtot_global;
 */

    my_malloc_init( rho_old, FP0_BASIS * ispin, REAL );

    my_malloc_init( sg_res, S0_BASIS, REAL );
#if USE_DIS_MAT
    sizeofmatrix = MXLLDA * MXLLDA;
#else
    sizeofmatrix = ct.num_states * ct.num_states;
#endif

#if !GAMMA_PT
    sizeofmatrix *= 2.0;
#endif

#if NONORTHO

    my_malloc_init( statearray, sizeofmatrix, REAL );
    my_malloc_init( l_s, sizeofmatrix, REAL );
    my_malloc_init( matB, sizeofmatrix, REAL );
    my_malloc_init( mat_X, sizeofmatrix, REAL );
    my_malloc_init( Hij, sizeofmatrix, REAL );
    my_malloc_init( theta, sizeofmatrix, REAL );
/*  my_malloc_init( work_dis, sizeofmatrix, REAL );
 * my_malloc_init( work_dis2, sizeofmatrix, REAL );
 * my_malloc_init( zz_dis, sizeofmatrix, REAL );
 * my_malloc_init( gamma_dis, sizeofmatrix, REAL );
 * my_malloc_init( uu_dis, sizeofmatrix, REAL );
*/

    my_malloc_init( work_matrix, ct.num_states * ct.num_states, REAL );

    /*  allocate memory for other uses  */

    item = ct.num_kpts * ct.num_states;
    item1 = max (ct.num_states * ct.num_states, item);
    item2 = 2 * P0_BASIS + S0_BASIS + item1;
    item = max (13 * S0_BASIS, item2);


    nproc = pct.nprow * pct.npcol;
    locr = ((ct.num_states / NB + 1) / nproc + 1) * NB + NB;
    lwork = 10 * (locr * 5 + NB);
    item1 = max (lwork, item);

    sl_init (&ictxt, pct.nprow, pct.npcol);
    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);
    if (myrow != -1)
        Cblacs_gridexit (ictxt);

    icrow = INDXG2P (&itwo, &nb, &myrow, &izero, &nprow);
    iccol = INDXG2P (&ione, &nb, &mycol, &izero, &npcol);
    mpc0 = NUMROC (&nn, &nb, &myrow, &icrow, &nprow);
    nqc0 = NUMROC (&nn, &nb, &mycol, &iccol, &npcol);
    nrc = NUMROC (&nn, &nb, &myrow, &izero, &npes);
    ldc = max (1, nrc);
    sizemqrleft = max ((NB * (NB - 1)) / 2, (nqc0 + mpc0) * NB) + NB * NB;
    sizemqrleft *= 2;
    qrmem = 2 * NN - 2;
    lwork = 5 * NN + NN * ldc + max (sizemqrleft, qrmem) + 1;
    item = max (lwork, item1);



#endif /* endif NONORTHO  */

}                               /* end allocate_matrix */

/********/
