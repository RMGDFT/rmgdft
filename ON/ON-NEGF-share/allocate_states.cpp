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
#include "prototypes_on.h"

void allocate_states()
{
    size_t size;
    size = (size_t)ct.num_states * (size_t)ct.num_kpts * (size_t)(ct.spin_flag +1);
    //printf("\n dsdsd %d %d %d", ct.num_states, ct.num_kpts, ct.spin_flag);
//    my_malloc(states, size, STATE);
    my_malloc(states1, size,  STATE);
    my_malloc(states_tem, size, STATE);
    my_malloc(state_to_proc, ct.num_states, int);

}

