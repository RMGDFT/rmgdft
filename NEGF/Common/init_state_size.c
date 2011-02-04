/************************** SVN Revision Information **************************
 **    $Id: init_state_size.c 1242 2011-02-02 18:55:23Z luw $    **
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
#include "md.h"

void init_state_size (STATE * states)
{

    int kpt, ispin, kpt1, st1;
    REAL *rptr, *rptr1, *rptr2;
    int state;
    int max_nx, max_ny, max_nz;


    max_nx = 0;
    max_ny = 0;
    max_nz = 0;

    if (pct.thispe == 0)
        printf ("\n  States orbital size: x.y.z   \n");

    for (state = 0; state < ct.num_states; state++)
    {
        states[state].orbit_nx = 2.0 * states[state].radius / ct.hxgrid / ct.xside + 3;
        states[state].orbit_ny = 2.0 * states[state].radius / ct.hygrid / ct.yside + 3;
        states[state].orbit_nz = 2.0 * states[state].radius / ct.hzgrid / ct.zside + 3;
        states[state].orbit_nx = min (states[state].orbit_nx, NX_GRID);
        states[state].orbit_ny = min (states[state].orbit_ny, NY_GRID);
        states[state].orbit_nz = min (states[state].orbit_nz, NZ_GRID);
        if (states[state].orbit_nx / 2 * 2 != states[state].orbit_nx)
            states[state].orbit_nx += 1;
        if (states[state].orbit_ny / 2 * 2 != states[state].orbit_ny)
            states[state].orbit_ny += 1;
        if (states[state].orbit_nz / 2 * 2 != states[state].orbit_nz)
            states[state].orbit_nz += 1;
        if (states[state].orbit_nx / 4 * 4 != states[state].orbit_nx)
            states[state].orbit_nx += 2;
        if (states[state].orbit_ny / 4 * 4 != states[state].orbit_ny)
            states[state].orbit_ny += 2;
        if (states[state].orbit_nz / 4 * 4 != states[state].orbit_nz)
            states[state].orbit_nz += 2;
        states[state].size = states[state].orbit_nx * states[state].orbit_ny
            * states[state].orbit_nz;
        max_nx = max (max_nx, states[state].orbit_nx);
        max_ny = max (max_ny, states[state].orbit_ny);
        max_nz = max (max_nz, states[state].orbit_nz);

#if 	LDEBUG
        if (pct.thispe == 0)
            printf (" %d: %d.%d.%d ", state, states[state].orbit_nx, states[state].orbit_ny,
                    states[state].orbit_nz);
        if (pct.thispe == 0 && (state - state / 5 * 5) == 0)
            printf ("\n");
#endif

    }

    ct.max_orbit_nx = max_nx;
    ct.max_orbit_ny = max_ny;
    ct.max_orbit_nz = max_nz;
    ct.max_orbit_size = max_nx * max_ny * max_nz;

    pct.psi_size = 0;
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        pct.psi_size += states[st1].size;
    }

    if (pct.thispe == 0)
        printf ("\n init_state_size Done! ");

}
