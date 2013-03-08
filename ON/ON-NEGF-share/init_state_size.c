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
#include "md.h"

void init_state_size(STATE * states)
{

    int kpt, ispin, kpt1, st1;
    REAL *rptr, *rptr1, *rptr2;
    int state;
    int max_nx, max_ny, max_nz;
    int nx_tem, ny_tem, nz_tem, item;

    max_nx = 0;
    max_ny = 0;
    max_nz = 0;

    if (pct.gridpe == 0)
        printf("\n  States orbital size: x.y.z   \n");

    item = 1<<ct.eig_parm.levels;
    for (state = 0; state < ct.num_states; state++)
    {
        nx_tem = 2.0 * states[state].radius / ct.hxgrid / ct.xside ;
        ny_tem = 2.0 * states[state].radius / ct.hygrid / ct.yside ;
        nz_tem = 2.0 * states[state].radius / ct.hzgrid / ct.zside ;
        //at coarsest level, the number of grid point in each direction
        //can be any number, at the upper levels, it must be an odd
        //number

        states[state].orbit_nx = (nx_tem + item -1)/item * item +1;
        states[state].orbit_ny = (ny_tem + item -1)/item * item +1;
        states[state].orbit_nz = (nz_tem + item -1)/item * item +1;


        states[state].size = states[state].orbit_nx * states[state].orbit_ny
            * states[state].orbit_nz;
        max_nx = max(max_nx, states[state].orbit_nx);
        max_ny = max(max_ny, states[state].orbit_ny);
        max_nz = max(max_nz, states[state].orbit_nz);
/*		if(pct.gridpe == 0) 
*			printf(" %d: %d.%d.%d ", 
*					state, states[state].orbit_nx, states[state].orbit_ny, 
*					states[state].orbit_nz); 
*		if( pct.gridpe == 0 && (state-state/5*5)==0 )
*			printf("\n"); 
*/
    }

    assert ( max_nx <= NX_GRID);
    assert ( max_ny <= NY_GRID);
    assert ( max_nz <= NZ_GRID);

    ct.max_orbit_nx = max_nx;
    ct.max_orbit_ny = max_ny;
    ct.max_orbit_nz = max_nz;
    ct.max_orbit_size = max_nx * max_ny * max_nz;

    pct.psi_size = 0;
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        pct.psi_size += states[st1].size;
    }

    if (pct.gridpe == 0)
        printf("\n init_state_size Done! ");

}
