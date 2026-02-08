#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
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
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "init_var.h"
#include "blacs.h"


void allocate_matrix_soft ()
{
    int ispin;
    

    int sbasis;

    ispin = ct.nspin;

    sbasis = (get_PX0_GRID() +2) * (get_PY0_GRID() +2) * (get_PZ0_GRID() +2);
    my_malloc_init( peaks, 100, double );

    my_malloc_init( rho, get_FP0_BASIS() * ispin, double );
    my_malloc_init( rhoc, get_FP0_BASIS(), double );
    my_malloc_init( vh, get_FP0_BASIS(), double );
    my_malloc_init( vnuc, get_FP0_BASIS(), double );
    my_malloc_init( vext, get_FP0_BASIS(), double );
    my_malloc_init( vcomp, get_FP0_BASIS(), double );
    my_malloc_init( vxc, get_FP0_BASIS() * ispin, double );
    my_malloc_init( vtot, get_FP0_BASIS(), double );
    my_malloc_init( vtot_c, get_P0_BASIS(), double );
    my_malloc_init( rhocore, get_FP0_BASIS(), double );
    
    if (ct.num_tfions > 0)
    {
        my_malloc_init( rho_tf, get_FP0_BASIS(), double );
    }

    my_malloc_init( rho_old, get_FP0_BASIS() * ispin, double );

    my_malloc_init( sg_res, sbasis, double );


}                               /* end allocate_matrix */

/********/
