/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/matrix_and_diag.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void matrix_and_diag(STATE *states)
 *   get the Hamiltonia and overlap matrixs and diagnolize it 
 * INPUTS
 *   states: a pointer for one kpoint state
 *   for each kpoint, we need do once.  
 * OUTPUT
 *   
 * PARENTS
 *   init.c scf.c
 * CHILDREN
 *   
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"


void matrix_and_diag(STATE * states, STATE * states1, rmg_double_t * vtot_c, int flag)
{

    int IA=1, JA=1, IB=1, JB=1, numst = ct.num_states;
    int level;




#if GAMMA_PT

    /* initialize matrices statearray and matB */
    get_HS(states, states1, vtot_c, Hij_00, Bij_00);

    Cpdgemr2d(numst, numst, Hij_00, IA, JA, pct.descb, Hij, IB, JB,
            pct.desca, pct.desca[1]);
    Cpdgemr2d(numst, numst, Bij_00, IA, JA, pct.descb, matB, IB, JB,
            pct.desca, pct.desca[1]);



    get_cholesky_real(matB);

    get_dm_diag_p(states, l_s, mat_X, Hij);
#else
    error_handler("not programmed for non-Gamma point");
#endif



}



