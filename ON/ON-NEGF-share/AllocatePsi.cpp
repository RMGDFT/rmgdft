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
#include <string.h>
#include <fcntl.h>
#include <assert.h>
#include <sys/mman.h>

#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "blas.h"
#include "RmgSumAll.h"

#include "prototypes_on.h"
#include "init_var.h"


void AllocatePsi(STATE * states, STATE * states1)
{

    int st1, idx ;
    size_t item;
    double *rptr, *rptr1, *rptr2, *rptr3;


    item = (ct.max_orbit_nx + 2) * (ct.max_orbit_ny + 2) * (ct.max_orbit_nz + 2);
    sg_orbit = new double[item];
    sg_orbit_res = new double[item];
    orbit_tem = new double[ct.max_orbit_size];

    size_t size, alloc_size;
    size = std::max(pct.psi_size, ct.state_per_proc * states[0].size);
    alloc_size = 3*size;
    rptr = new double[alloc_size];
    if(NULL == rptr) 
    {
        printf("\n cannot locate memory for psi %lu %lu \n", size, alloc_size);
        exit(0);
    }

    rptr1 = rptr + size;
    rptr2 = rptr + 2*size;

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

    
     if(ct.LocalizedOrbitalLayout == LO_projection)
         return;
//  count how many orbtials will be received from other processor and 
//    allocate memory for them in states[].psiR

    int state_per_proc = ct.state_per_proc + 2;

    int i, loop, num_recv;
    size_t tot_recv;
    
    tot_recv = 0;
    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {
        num_recv = recv_from1[loop * state_per_proc + 1];
        tot_recv += num_recv;
    }


    ct.num_orbitals_total = ct.state_end - ct.state_begin + tot_recv;

    ct.orbitals_list = new int[ct.num_orbitals_total];

    idx = 0;
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        ct.orbitals_list[idx] = st1;
        idx++;
    }

    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {
        num_recv = recv_from1[loop * state_per_proc + 1];
        for(i = 0; i< num_recv; i++)
        {
            st1 = recv_from1[loop * state_per_proc + i + 2];
            ct.orbitals_list[idx] = st1;
            idx++;

        }

    }

    
    item = ct.max_orbit_nx  * ct.max_orbit_ny  * ct.max_orbit_nz;

    alloc_size =  (size_t)tot_recv * item ;

    ct.nvme_orbital_fd = -1;
    // orbitals from other proces  are actually stored here

    if(ct.nvme_orbitals)
    {
        std::string newpath;
        if(ct.nvme_orbital_fd != -1) close(ct.nvme_orbital_fd);

        newpath = ct.nvme_orbitals_path + std::string("rmg_orbital") + std::to_string(pct.spinpe) +
            std::to_string(pct.kstart) + std::to_string(pct.gridpe);
        ct.nvme_orbital_fd = FileOpenAndCreate(newpath, O_RDWR|O_CREAT|O_TRUNC, (mode_t)0600);
        rptr3 = (double *)CreateMmapArray(ct.nvme_orbital_fd, alloc_size * sizeof(double));
        if(!rptr3) rmg_error_handler(__FILE__,__LINE__,"Error: CreateMmapArray failed for orbitals. \n");
        madvise(rptr3, alloc_size * sizeof(double), MADV_SEQUENTIAL);
    }
    else
    {
        rptr3 = new double[alloc_size];
    }



    if(NULL == rptr3) 
    {
        printf("\n cannot locate memory for recv %lu %lu %lu \n", item, tot_recv, alloc_size);
        exit(0);
    }

    //dprintf("\n tot_recv  %d", tot_recv);


    // only assign memory for orbitals on other procs.
    for (st1 = 0; st1 < ct.num_orbitals_total; st1++)
    {
        idx = ct.orbitals_list[st1];
        if(idx < ct.state_begin || idx >= ct.state_end)
        {
            states[idx].psiR = rptr3;
            rptr3 += item;
        }
    }

}

