/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/allocate_psi.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void 
 *        
 * INPUTS
 *   
 * OUTPUT
 *
 * PARENTS
 *   init.c
 * CHILDREN
 * 
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"

void allocate_psi(STATE * states, STATE * states1)
{

    int kpt, ispin, kpt1, st1, item;
    rmg_double_t *rptr, *rptr1, *rptr2;


    item = (ct.max_orbit_nx + 2) * (ct.max_orbit_ny + 2) * (ct.max_orbit_nz + 2);
    my_malloc_init( sg_orbit, item, rmg_double_t );
    my_malloc_init( sg_orbit_res, item, rmg_double_t );
    my_malloc_init( orbit_tem, ct.max_orbit_size, rmg_double_t );

    my_malloc_init( rptr, pct.psi_size, rmg_double_t );
    my_malloc_init( rptr1, pct.psi_size, rmg_double_t );
    my_malloc_init( rptr2, pct.psi_size, rmg_double_t );

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        states[st1].psiR = rptr;
        states1[st1].psiR = rptr1;
        states_tem[st1].psiR = rptr2;
        rptr += states[st1].size;
        rptr1 += states[st1].size;
        rptr2 += states[st1].size;
    }

    if (pct.gridpe == 0)
        printf("\n allocate_psi Done! ");

}
