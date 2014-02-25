/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/init_soft.c *****
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
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"

#define min(a,b) (((a)>(b)) ? (b) : (a))



void init_soft (rmg_double_t * vh, rmg_double_t * rho, rmg_double_t * rhocore, rmg_double_t * rhoc,
                STATE * states, STATE * states1, rmg_double_t * vnuc, rmg_double_t * vext, rmg_double_t * vxc, rmg_double_t * vh_old,
                rmg_double_t * vxc_old, STATE *states_distribute)
{

    int kpt, ic, idx, ion, ispin, kpt1;
    int flag, level;
    rmg_double_t time1, time2, cut_init;
    rmg_double_t tem;
    int item;
    char *nameL, *nameC, *nameR;
    int st1, iprobe, i;
        int ione = 1;


    time1 = my_crtc ();

    ct.psi_nbasis = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();
    ct.psi_nxgrid = get_NX_GRID();
    ct.psi_nygrid = get_NY_GRID();
    ct.psi_nzgrid = get_NZ_GRID();

    ct.psi_fnbasis = get_FNX_GRID() * get_FNY_GRID() * get_FNZ_GRID();
    ct.psi_fnxgrid = get_FNX_GRID();
    ct.psi_fnygrid = get_FNY_GRID();
    ct.psi_fnzgrid = get_FNZ_GRID();



    flag = 0;


    ct.states = states;

    /* initialize the lattice basis vectors */
    int ibrav = get_ibrav_type();
    latgen (&ibrav, ct.celldm, ct.a0, ct.a1, ct.a2, &ct.omega, &flag);

    init_parameter (states);

    /* initialize the reciprocal lattice vectors */
    recips ();


    /* Get the crystal or cartesian coordinates of the ions */
    init_pos ();


    /* Set initial ionic coordinates to the current ones. */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        for (ic = 0; ic < 3; ic++)
        {
            ct.ions[ion].icrds[ic] = ct.ions[ion].crds[ic];
            ct.ions[ion].ixtal[ic] = ct.ions[ion].xtal[ic];
        }                       /* end for ic */

    }                           /* end for */



    if (pct.gridpe == 0)
    {

        /* Write header to stdout */
        write_header ();

    }                           /* endif */

    my_barrier ();

    /* Initialize the mehrstellen weights */
    /*get_mehr (); */

    get_state_to_proc (states);

    init_state_size (states);

    state_corner_xyz (states);


    /* allocate memory for wave functions states.psiR and psiI */
    allocate_psi (states, states1);
    if (pct.gridpe == 0)
        printf ("Allocate_psi is done \n");
    fflush (NULL);


    is_state_overlap (states, state_overlap_or_not);
    get_orbit_overlap_region (states);
    init_comm (states);

    duplicate_states_info (states, states1);
    my_barrier ();


    pmo_init();
/*
    nameL = lcr[1].name;
    nameC = lcr[0].name;
    nameR = lcr[2].name;
*/
        read_orbital(states);
    interpolation_orbit (states);
    init_state_distribute(states, states_distribute);

    scale_orbital(states, states_distribute);
#if GPU_ENABLED
    init_gpu();
    cublasSetVector( pct.num_local_orbit * get_P0_BASIS(), sizeof( double ), states_distribute[0].psiR, ione, ct.gpu_states, ione );

#endif
        if (pct.gridpe == 0) printf ("completed: read_orbital \n");
        allocate_matrix_LCR();
        if (pct.gridpe == 0) printf ("completed: allocate_matrix \n");
    if (ct.runflag > 111)
    {
        read_lead_matrix();
        if (pct.gridpe == 0) printf ("completed: read_lead_matrix \n");
    }

    /*exit(0); */ 

/*
    nameL = lcr[1].lead_name;
    nameC = lcr[0].lead_name;
    nameR = lcr[2].lead_name;
*/
    if(ct.runflag <113)
    {
        read_potrho_LCR (vh, vxc, rho);
        if (pct.gridpe == 0) printf ("completed: read_potrho_LCR \n");
    }
    else
    {    
        read_rho_and_pot (ct.infile, vh, vxc, vh_old, vxc_old, rho);
        if (pct.gridpe == 0) printf ("completed: read_rho_and_pot \n");
    }
     

    if(ct.runflag == 300) 
    {
        plane_average_rho(rho); 
        exit(0);
    }

    write_rho_x (vh, "vh_init");  // information about vh_init is recorded in zvec array

    if (ct.vcomp_Rend > 0) // if ct.vcomp_Rend > 0 means the user turned on the initial compensating potential correction
    {
	    init_comp (vh);
	    for (idx = 0; idx < get_FP0_BASIS(); idx++)
	    {
		    vh[idx] = vh[idx] + vcomp[idx]; //add compensating potential to align lead and center part at the very beginning
	    }

	    write_rho_x (vh, "vh_vcomp_init");
    }

/*  interpolation for the orbits if the grid space is slightly different between lead and conductor */


    allocate_masks (states);

    for (level = 0; level < ct.eig_parm.levels + 1; level++)
        make_mask_grid_state (level, states);


/*    normalize_orbits(states);
*/
//    ortho_norm_local (states);


    /*modify the lead Hamiltonia by bias and Fermi energy */
    for(iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        tem = (-lcr[iprobe].EF_old + lcr[iprobe].bias + lcr[iprobe].EF_new) * eV_Ha;

        idx = pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];

        for (st1 = 0; st1 < idx; st1++)
        {
            lcr[iprobe].H00[st1] += lcr[iprobe].S00[st1] * tem;
            lcr[iprobe].H01[st1] += lcr[iprobe].S01[st1] * tem;
        }


		i = cei.probe_in_block[iprobe - 1];
        idx = pmo.mxllda_cond[i] * pmo.mxlocc_lead[iprobe-1];

        for (st1 = 0; st1 < idx; st1++)
        {
            lcr[iprobe].HCL[st1] += lcr[iprobe].SCL[st1] * tem;
        }


    }
    my_barrier ();


    for (level = 0; level < ct.eig_parm.levels + 1; level++)
        make_mask_grid_state (level, states);

    /* Initialize the radial potential stuff */
    init_psp_soft ();

    /* Initialize the radial qfunction stuff */
    init_qfunct ();

    /* Initialize symmetry stuff */
    //init_sym ();

    /* Initialize the nuclear local potential and the compensating charges */
    init_nuc (vnuc, rhoc, rhocore);

    init_ext (vext, ct.gbias_begin, ct.gbias_end, ct.BT, ct.gate_bias);
 
    write_rho_x (rho, "rho_init1");


    /* Initialize Non-local operators */
    init_nl_xyz ();
    get_ion_orbit_overlap_nl (states);

    get_nlop ();
    init_nonlocal_comm ();

    /* Initialize qfuction in Cartesin coordinates */
    get_QI ();

    /* Get the qqq */
    get_qqq ();


    time2 = my_crtc ();
    rmg_timings (INIT_SOFT_TIME, time2 - time1);



}                               /* end init */

