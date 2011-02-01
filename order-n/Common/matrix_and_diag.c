/************************** SVN Revision Information **************************
 **    $Id: matrix_and_diag.c 587 2006-08-18 22:57:56Z miro $    **
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
#include "md.h"



void matrix_and_diag(STATE * states, STATE * states1, REAL * vtot_c, int flag)
{

    char side = 'l', uplo = 'l';
    double zero = 0.0, one = 1.0;
    int numst = ct.num_states;
    int maxst = numst;
    int level;
    REAL tem, tem1, tem2;
    int idx, idx1;
    double time1, time2;



#if GAMMA_PT

    /* initialize matrices statearray and matB */
    time1 = my_crtc();
    get_HS(states, states1, vtot_c, Hij, matB);

    time2 = my_crtc();
    rmg_timings(HIJ_TIME, time2 - time1, 0);


    get_cholesky_real(matB);

    get_dm_diag_p(states, l_s, mat_X, Hij);
#else
    error_handler("not programmed for non-Gamma point");
#endif



}
