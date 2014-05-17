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
#include "init_var.h"


void allocate_states()
{
    int size;
    my_malloc(states, ct.num_states, STATE);
    my_malloc(states1, ct.num_states, STATE);
    my_malloc(states_tem, ct.num_states, STATE);
    my_malloc(states_distribute, ct.num_states, STATE);
    my_malloc(state_to_ion, ct.num_states, int);
    my_malloc(state_to_proc, ct.num_states, int);

    my_malloc(vloc_state_overlap_or_not, ct.num_ions * ct.num_states, char);
}

