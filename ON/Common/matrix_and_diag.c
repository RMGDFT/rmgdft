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


void matrix_and_diag(STATE * states, STATE * states1, double * vtot_c, int flag)
{

    int IA=1, JA=1, IB=1, JB=1, numst = ct.num_states;
    int level;




#if GAMMA_PT

    /* initialize matrices statearray and matB */
void *RT = BeginRmgTimer("3-matrix_and_diag");

void *RT1 = BeginRmgTimer("3-matrix_and_diag: get_HS");
    get_HS(states, states1, vtot_c, Hij_00, Bij_00);
    EndRmgTimer(RT1);
//    DiagElemental(numst, Hij_00, Bij_00);

void *RT2 = BeginRmgTimer("3-matrix_and_diag: cpdgemr2d");
    Cpdgemr2d(numst, numst, Hij_00, IA, JA, pct.descb, Hij, IB, JB,
           pct.desca, pct.desca[1]);
   Cpdgemr2d(numst, numst, Bij_00, IA, JA, pct.descb, matB, IB, JB,
          pct.desca, pct.desca[1]);

    EndRmgTimer(RT2);


//    get_cholesky_real(matB);

void *RT3 = BeginRmgTimer("3-matrix_and_diag: diag");
    get_dm_diag_p(states, matB, mat_X, Hij);
    EndRmgTimer(RT3);

    EndRmgTimer(RT);

#else
    error_handler("not programmed for non-Gamma point");
#endif



}



