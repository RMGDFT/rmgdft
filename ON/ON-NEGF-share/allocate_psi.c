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
#include "prototypes_on.h"
#include "init_var.h"

void allocate_psi(STATE * states, STATE * states1)
{

    int st1, item;
    double *rptr, *rptr1, *rptr2, *rptr3;


    item = (ct.max_orbit_nx + 2) * (ct.max_orbit_ny + 2) * (ct.max_orbit_nz + 2);
    my_malloc_init( sg_orbit, item, double );
    my_malloc_init( sg_orbit_res, item, double );
    my_malloc_init( orbit_tem, ct.max_orbit_size, double );

    my_malloc_init( rptr, pct.psi_size, double );
    my_malloc_init( rptr1, pct.psi_size, double );
    my_malloc_init( rptr2, pct.psi_size, double );

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

//  count how many orbtials will be received from other processor and 
//    allocate memory for them in states[].psiR

    int state_per_proc = ct.state_per_proc + 2;

    int i, loop, num_recv, tot_recv;
    
    tot_recv = 0;
    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {
        num_recv = recv_from1[loop * state_per_proc + 1];
        tot_recv += num_recv;
    }

    dprintf("\n total recv orbitals %d", tot_recv);
    item = ct.max_orbit_nx  * ct.max_orbit_ny  * ct.max_orbit_nz;
    my_malloc_init(rptr3, tot_recv * item, double);
    dprintf("\n tot_recv  %d", tot_recv);

    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {
        num_recv = recv_from1[loop * state_per_proc + 1];
        for(i = 0; i< num_recv; i++)
        {
            st1 = recv_from1[loop * state_per_proc + i + 2];
            states[st1].psiR = rptr3;
            rptr3  += item;
        }

    }


}
