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
 *   void init(REAL *vh, REAL *rho, REAL *rhocore, REAL *rhoc, 
 *             STATE *states, REAL *vnuc, REAL *vxc)
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
#include "md.h"

#define min(a,b) (((a)>(b)) ? (b) : (a))

void init_wf_atom(STATE *);

void init(REAL * vh, REAL * rho, REAL * rhocore, REAL * rhoc,
          STATE * states, STATE * states1, REAL * vnuc, REAL * vxc, REAL * vh_old, REAL * vxc_old)
{

    int kpt, ic, idx, ion, ispin, kpt1;
    int flag, level;
    REAL time1, time2, cut_init;
    double tem, tem1;

   int  gridpe = pct.gridpe;

    time1 = my_crtc();
    /* initialize the lattice basis vectors */
    flag = 0;

    ct.psi_nbasis = NX_GRID * NY_GRID * NZ_GRID;
    ct.psi_nxgrid = NX_GRID;
    ct.psi_nygrid = NY_GRID;
    ct.psi_nzgrid = NZ_GRID;

    ct.psi_fnbasis = FNX_GRID * FNY_GRID * FNZ_GRID;
    ct.psi_fnxgrid = FNX_GRID;
    ct.psi_fnygrid = FNY_GRID;
    ct.psi_fnzgrid = FNZ_GRID;



    my_malloc_init( ct.energies, ct.max_scf_steps, REAL );

    ct.states = states;

    init_parameter(states);
    if(gridpe == 0) printf("\n init_parameter done %f sec",my_crtc()-time1 );

    latgen(&ct.ibrav, ct.celldm, ct.a0, ct.a1, ct.a2, &ct.omega, &flag);
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

    get_state_to_proc(states);

    if(gridpe == 0) printf("\n init get_state_to_proc done %f sec",my_crtc()-time1 );

    /* allocate memory for wave functions states.psiR and psiI */
    init_state_size(states);
    if(gridpe == 0) printf("\n init_state_size  done %f sec",my_crtc()-time1 );

    state_corner_xyz(states);
    if(gridpe == 0) printf("\n init_state_corner_xyz  done %f sec",my_crtc()-time1 );

    allocate_psi(states, states1);
    if(gridpe == 0) printf("\n init_allocate_psi  done %f sec",my_crtc()-time1 );

    init_states();

    is_state_overlap(states);

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

    if (ct.runflag == 0)
    {
        cut_init = 1.0;
        make_mask_grid(cut_init, 0, states);

        /* Set initial states to random start */
        for (ispin = 0; ispin <= ct.spin; ispin++)
            for (kpt = pct.kstart; kpt < pct.kend; kpt++)
            {
                kpt1 = kpt + ispin * ct.num_kpts;
                init_wf(&states[kpt1 * ct.num_states]);
            }

    if(gridpe == 0) printf("\n init_wf done %f sec",my_crtc()-time1 );
    }

    if (ct.runflag == INIT_GAUSSIAN)
    {

        /* Set initial states to random start */
        for (ispin = 0; ispin <= ct.spin; ispin++)
            for (kpt = pct.kstart; kpt < pct.kend; kpt++)
            {
                kpt1 = kpt + ispin * ct.num_kpts;
                init_wf_gaussian(&states[kpt1 * ct.num_states]);
            }

    if(gridpe == 0) printf("\n init_wf_gaussian done %f sec",my_crtc()-time1 );
    }


    for (level = 0; level < ct.eig_parm.levels + 1; level++)
        make_mask_grid_state(level, states);

    /* Initialize the radial potential stuff */
    init_kbr();
    if(gridpe == 0) printf("\n init_kbr done %f sec",my_crtc()-time1 );

    /* Initialize symmetry stuff */
    init_sym();
    if(gridpe == 0) printf("\n init_sys done %f sec",my_crtc()-time1 );

    /* Initialize the nuclear local potential and the compensating charges */
    init_nuc(vnuc, rhoc, rhocore);
    if(gridpe == 0) printf("\n init_nuc done %f sec",my_crtc()-time1 );

    if (ct.runflag == INIT_FIREBALL)
    {
        init_wf_atom(states);
        for (idx = 0; idx < FP0_BASIS; idx++)
            vh[idx] = ZERO;
        for (idx = 0; idx < FP0_BASIS; idx++)
            vxc[idx] = ZERO;

        init_rho_atom(rho);
    if(gridpe == 0) printf("\n init_rho_atom done %f sec",my_crtc()-time1 );

#if DEBUG
        tem = 0.0;
        tem1 = 0.0;

        for (idx = 0; idx < FP0_BASIS; idx++)

        {
            tem += rho[idx];
            tem1 += rhoc[idx];
        }
        tem = real_sum_all(tem, pct.grid_comm);
        tem1 = real_sum_all(tem1, pct.grid_comm);
        if (pct.gridpe == 0)
            printf("\n %f %f initial rho sum ", tem, tem1);
        write_rho_x(rho, "rhoooo");
        write_rho_x(rhoc, "rhoccc");
#endif

        get_vxc(rho, rhocore, vxc);
        pack_vhstod(vh, ct.vh_ext, FPX0_GRID, FPY0_GRID, FPZ0_GRID);
        get_vh(rho, rhoc, vh, 30, ct.poi_parm.levels);
        for (idx = 0; idx < FP0_BASIS; idx++)
            vh_old[idx] = vh[idx];
        for (idx = 0; idx < FP0_BASIS; idx++)
            vxc_old[idx] = vxc[idx];
    if(gridpe == 0) printf("\n init_vpot done %f sec",my_crtc()-time1 );

    }

    if (ct.runflag == 0 | ct.runflag == INIT_GAUSSIAN)
    {
        /* Set the initial hartree potential to a constant */
        for (idx = 0; idx < FP0_BASIS; idx++)
            vh[idx] = ZERO;
        for (idx = 0; idx < FP0_BASIS; idx++)
            vxc[idx] = ZERO;
        for (idx = 0; idx < FP0_BASIS; idx++)
            rho[idx] = rhoc[idx];

        get_vxc(rho, rhocore, vxc);
        pack_vhstod(vh, ct.vh_ext, FPX0_GRID, FPY0_GRID, FPZ0_GRID);
        get_vh(rho, rhoc, vh, 3, ct.poi_parm.levels);
        for (idx = 0; idx < FP0_BASIS; idx++)
            vh_old[idx] = vh[idx];
        for (idx = 0; idx < FP0_BASIS; idx++)
            vxc_old[idx] = vxc[idx];
    if(gridpe == 0) printf("\n init_vpot_0 done %f sec",my_crtc()-time1 );

        if (ct.spin == 1)
        {
            for (idx = 0; idx < FP0_BASIS; idx++)
            {
                rho[idx + FP0_BASIS] = 0.4 * rho[idx];
                rho[idx] = 0.6 * rho[idx];
            }
        }
    }


    /* Initialize Non-local operators */
    init_nl_xyz();
    if(gridpe == 0) printf("\n init_nl_xyz done %f sec",my_crtc()-time1 );
    get_ion_orbit_overlap_nl(states);
    if(gridpe == 0) printf("\n init_get_ion_orbit_overlap_nl done %f sec",my_crtc()-time1 );

    get_nlop();
    if(pct.gridpe == 0) printf("\n pe %d init_get_lop done %f sec",pct.gridpe, my_crtc()-time1 );

    /* If not an initial run read data from files */
    if (ct.runflag == 1)
    {
        read_data(ct.infile, vh, vxc, vh_old, vxc_old, rho, states);
    if(gridpe == 0) printf("\n init_read_data done %f sec",my_crtc()-time1 );
        pack_vhstod(vh, ct.vh_ext, FPX0_GRID, FPY0_GRID, FPZ0_GRID);
    }

    my_barrier();
    if(pct.gridpe == 0) printf("\n pe %d my_barrier %f sec",pct.gridpe, my_crtc()-time1 );

    init_nonlocal_comm();
    if(pct.gridpe == 0) printf("\n pe %d init_nonlocal_comm done %f sec",pct.gridpe,my_crtc()-time1 );

    /* Initialize qfuction in Cartesin coordinates */
    init_qfunct();
    if(gridpe == 0) printf("\n init_qfunct done %f sec",my_crtc()-time1 );
    get_QI();
    if(gridpe == 0) printf("\n init_get_QI done %f sec",my_crtc()-time1 );

    /* Get the qqq */
    get_qqq();
    if(gridpe == 0) printf("\n init_get_qqq done %f sec",my_crtc()-time1 );


    /* If diagonalization is requested do a subspace diagonalization */
    for (ispin = 0; ispin <= ct.spin; ispin++)
    {
        for (kpt = pct.kstart; kpt < pct.kend; kpt++)
        {
            kpt1 = kpt + ispin * ct.num_kpts;
            for (idx = 0; idx < FP0_BASIS; idx++)
                vtot[idx] = vh_old[idx] + vxc_old[idx] + vnuc[idx];

            get_vtot_psi(vtot_c, vtot, FG_NX);
            get_ddd(vtot);

            flag = 0;
            matrix_and_diag(ct.kp[kpt1].kstate, states1, vtot_c, flag);
    if(gridpe == 0) printf("\n init_matrix_and_diag done %f sec",my_crtc()-time1 );
        }
        if (pct.gridpe == 0)
            write_eigs(states);
    if(gridpe == 0) printf("\n init_write_eigs done %f sec",my_crtc()-time1 );
    }

    time2 = my_crtc();
    rmg_timings(INIT_TIME, (time2 - time1));
    my_barrier();
    if(gridpe == 0) printf("\n init_barrier done %f sec",my_crtc()-time1 );


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
