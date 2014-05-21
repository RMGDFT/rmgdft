/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/init.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void init(rmg_double_t *vh, rmg_double_t *rho, rmg_double_t *rhocore, rmg_double_t *rhoc, 
 *             STATE *states, rmg_double_t *vnuc, rmg_double_t *vxc)
 *   Basic initialization stuff.
 * INPUTS
 *   rhocore: core charge density for non-linear core corection
 *            it was read in md.c by call read_pseudo.c
 * OUTPUT
 *   vh:  initial Hartree potential 
 *   rho: initial charge density
 *   rhoc: compensating charge density
 *   states: wave functions are initialize randomly or read from data
 *   vnuc:  pseudopotential
 *   vxc: exchange correlation potential calculated from initial rho
 * PARENTS
 *   run.c
 * CHILDREN
 *   latgen.c recips.c init_pos.c read_data.c init_wf.c init_kbr.c
 *   init_sys.c get_nlop.c init_nuc.c init_wflcao.c get_vxc.c
 *   pack_vhstod.c 
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
#define min(a,b) (((a)>(b)) ? (b) : (a))

void init_wf_atom(STATE *);

void init(rmg_double_t * vh, rmg_double_t * rho, rmg_double_t * rhocore, rmg_double_t * rhoc,
          STATE * states, STATE * states1, rmg_double_t * vnuc, rmg_double_t * vxc, rmg_double_t * vh_old, rmg_double_t * vxc_old)
{

    int ic, idx, ion;
    int level;
    rmg_double_t time1, time2;

   int  gridpe = pct.gridpe;

    time1 = my_crtc();
    /* initialize the lattice basis vectors */

    ct.psi_nbasis = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();

    ct.psi_fnbasis = get_FNX_GRID() * get_FNY_GRID() * get_FNZ_GRID();


    my_malloc_init( ct.energies, ct.max_scf_steps, rmg_double_t );

    ct.states = states;

    init_parameter(states);
    if(gridpe == 0) printf("\n init_parameter done %f sec",my_crtc()-time1 );

    if(gridpe == 0) printf("\n init latgen done %f sec",my_crtc()-time1 );

    init_parameter(states);
    if(gridpe == 0) printf("\n init_parameter done %f sec",my_crtc()-time1 );

    /* initialize the reciprocal lattice vectors */
    recips();

    /* Get the crystal or cartesian coordinates of the ions */
    init_pos();

    /* Set initial ionic coordinates to the current ones. */
    for (ion = 0; ion < ct.num_ions; ion++)
        for (ic = 0; ic < 3; ic++)
        {
            ct.ions[ion].icrds[ic] = ct.ions[ion].crds[ic];
            ct.ions[ion].ixtal[ic] = ct.ions[ion].xtal[ic];
        }

    if (pct.gridpe == 0)
    {
        write_header();
    }

    if(gridpe == 0) printf("\n init write_header done %f sec",my_crtc()-time1 );
    my_barrier();
    if(gridpe == 0) printf("\n init my_brarrier done %f sec",my_crtc()-time1 );

    /* Initialize the mehrstellen weights */
    get_mehr();

    if(gridpe == 0) printf("\n init get_state_to_proc done %f sec",my_crtc()-time1 );

    /* allocate memory for wave functions states.psiR and psiI */
    init_state_size(states);
    if(gridpe == 0) printf("\n init_state_size  done %f sec",my_crtc()-time1 );

    state_corner_xyz(states);
    if(gridpe == 0) printf("\n init_state_corner_xyz  done %f sec",my_crtc()-time1 );

    allocate_psi(states, states1);


    if(gridpe == 0) printf("\n init_allocate_psi  done %f sec",my_crtc()-time1 );

    init_states();

    is_state_overlap(states, state_overlap_or_not);

    if(gridpe == 0) printf("\n init_is_state_overlap  done %f sec",my_crtc()-time1 );

    get_orbit_overlap_region(states);
    if(gridpe == 0) printf("\n init_get_orbit_overlap_region  done %f sec",my_crtc()-time1 );

    init_comm(states);

    init_comm_res(states);
    if(gridpe == 0) printf("\n init_comm_ress  done %f sec",my_crtc()-time1 );

    duplicate_states_info(states, states1);
    duplicate_states_info(states, states_tem);

    if(gridpe == 0) printf("\n init_duplicate_states_info done %f sec",my_crtc()-time1 );
    my_barrier();
    if(gridpe == 0) printf("\n init_barrier done %f sec",my_crtc()-time1 );


    allocate_masks(states);
    if(gridpe == 0) printf("\n init_allocate mask done %f sec",my_crtc()-time1 );

    for (level = 0; level < ct.eig_parm.levels + 1; level++)
        make_mask_grid_state(level, states);

    /* Initialize the radial potential stuff */
    init_kbr();
    if(gridpe == 0) printf("\n init_kbr done %f sec",my_crtc()-time1 );

    /* Initialize symmetry stuff */
   //  init_sym();
    if(gridpe == 0) printf("\n init_sys done %f sec",my_crtc()-time1 );

    /* Initialize the nuclear local potential and the compensating charges */
    init_nuc(vnuc, rhoc, rhocore);
    if(gridpe == 0) printf("\n init_nuc done %f sec",my_crtc()-time1 );



    /* Initialize Non-local operators */
    init_nl_xyz();
    if(gridpe == 0) printf("\n init_nl_xyz done %f sec",my_crtc()-time1 );
    get_ion_orbit_overlap_nl(states);
    if(gridpe == 0) printf("\n init_get_ion_orbit_overlap_nl done %f sec",my_crtc()-time1 );

    get_nlop();
    if(pct.gridpe == 0) printf("\n pe %d init_get_lop done %f sec",pct.gridpe, my_crtc()-time1 );

    fflush(NULL);

    my_barrier();
    if(pct.gridpe == 0) printf("\n pe %d my_barrier %f sec",pct.gridpe, my_crtc()-time1 );
    fflush(NULL);

    init_nonlocal_comm();
    if(pct.gridpe == 0) printf("\n pe %d init_nonlocal_comm done %f sec",pct.gridpe,my_crtc()-time1 );
    fflush(NULL);

    /* Initialize qfuction in Cartesin coordinates */
    init_qfunct();
    if(gridpe == 0) printf("\n init_qfunct done %f sec",my_crtc()-time1 );
    fflush(NULL);
    get_QI();
    if(gridpe == 0) printf("\n init_get_QI done %f sec",my_crtc()-time1 );
    fflush(NULL);

    /* Get the qqq */
    get_qqq();
    if(gridpe == 0) printf("\n init_get_qqq done %f sec",my_crtc()-time1 );
    fflush(NULL);

    for (idx = 0; idx < get_FP0_BASIS(); idx++) vh[idx] = ZERO;
    switch(ct.runflag)
    {
        case 0:
            init_wf(states);

            for (idx = 0; idx < get_FP0_BASIS(); idx++)
                rho[idx] = rhoc[idx];
            break;
        case INIT_FIREBALL:
            init_wf_atom(states);
            init_rho_atom(rho);
            printf("\n init_rho_atom done %f sec",my_crtc()-time1 );
            break;

        case INIT_GAUSSIAN:
            init_wf_gaussian(states);
            for (idx = 0; idx < get_FP0_BASIS(); idx++)
                rho[idx] = rhoc[idx];
            break;

        case 1:
            time1 = my_crtc();
            read_data(ct.infile, vh, vxc, vh_old, vxc_old, rho, states);
            if(gridpe == 0) printf("\n init_read_data done %f sec",my_crtc()-time1 );
            pack_vhstod(vh, ct.vh_ext, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), ct.boundaryflag);
            break;

    }

    if(ct.runflag !=1) 
    {
        get_vxc(rho, rho, rhocore, vxc);
        pack_vhstod(vh, ct.vh_ext, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), ct.boundaryflag);
        get_vh (rho, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, 0.0, ct.boundaryflag);
        for (idx = 0; idx < get_FP0_BASIS(); idx++)
            vh_old[idx] = vh[idx];
        for (idx = 0; idx < get_FP0_BASIS(); idx++)
            vxc_old[idx] = vxc[idx];
    }


    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx];

    get_vtot_psi(vtot_c, vtot, get_FG_RATIO());
    get_ddd(vtot);


    time2 = my_crtc();
    rmg_timings(INIT_TIME, (time2 - time1));
    my_barrier();
    fflush(NULL);


#if	DEBUG
    print_state_sum(states);
    print_status(states, vh, vxc, vnuc, vh_old, "before leaving init.c  ");
    print_sum(pct.psi_size, states[ct.state_begin].psiR, "init.c states sum ");
#endif

    /* some utilities, used in debuging */

    if (pct.gridpe > 100000)
    {
        print_status(states, vh, vxc, vnuc, vh_old, "before leaving init.c  ");
        print_state_sum(states1);
        print_state(&states[0]);
        print_sum(pct.psi_size, states[ct.state_begin].psiR, "init.c states sum ");
    }


}
