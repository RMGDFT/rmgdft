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
 *   void init(double *vh, double *rho, double *rhocore, double *rhoc, 
 *             STATE *states, double *vnuc, double *vxc)
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
 *   init_sys.c init_nuc.c init_wflcao.c get_vxc.c
 *   pack_vhstod.c 
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "LCR.h"
#include "prototypes_on.h"
#include "prototypes_negf.h"
#include "init_var.h"



#include "Scalapack.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"

#include "pmo.h"

#include "RmgException.h"
#include "RmgParallelFft.h"
#include "Atomic.h"



void is_state_overlap (STATE *, char *);

void InitNegf (double * vh, double * rho, double * rhocore, double * rhoc, double * rho_tf,
                STATE * states, STATE * states1, double * vnuc, double * vext, double * vxc, double * vh_old,
                double * vxc_old, std::unordered_map<std::string, InputKey *>& ControlMap)
{

    int kpt, ic, idx, ion, ispin, kpt1;
    int flag, level;
    double cut_init;
    double tem;
    int item;
    char *nameL, *nameC, *nameR;
    int st1, iprobe, i;
        int ione = 1;



    ct.psi_nbasis = Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1);



    flag = 0;


    ct.states = states;

    /* initialize the lattice basis vectors */

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


    if (ct.num_tfions > 0)
    {
	get_tf_rho(rho_tf);
    }



    if (pct.gridpe == 0)
    {

        /* Write header to stdout */
        write_header ();

    }                           /* endif */

    MPI_Barrier(pct.img_comm);

    /* Initialize the mehrstellen weights */
    /*get_mehr (); */


    init_state_size (states);

    state_corner_xyz (states);



    int size = (ct.state_end - ct.state_begin) * ct.num_states;

    state_overlap_or_not = new char[size];


    is_state_overlap (states, state_overlap_or_not);
    get_orbit_overlap_region (states);
    GetOrbitalPairs(states);
    init_comm (states);
    init_comm_res (states);

    AllocatePsi (states, states1);

    duplicate_states_info (states, states1);
    MPI_Barrier(pct.img_comm);


    pmo_init();
/*
    nameL = lcr[1].name;
    nameC = lcr[0].name;
    nameR = lcr[2].name;
*/

    RmgTimer *RT1 = new RmgTimer("1-TOTAL: init:  read_orbital");
    read_orbital(states);
    interpolation_orbit (states);
    delete(RT1);

    scale_orbital(states);

    if (pct.imgpe == 0) printf ("completed: read_orbital \n");
    allocate_matrix_LCR();
    if (pct.imgpe == 0) printf ("completed: allocate_matrix \n");
    if (ct.runflag > 111)
    {
        read_lead_matrix();
        if (pct.imgpe == 0) printf ("completed: read_lead_matrix \n");
    }

    /*exit(0); */ 

    /*
       nameL = lcr[1].lead_name;
       nameC = lcr[0].lead_name;
       nameR = lcr[2].lead_name;
     */
    RmgTimer *RT3 = new RmgTimer("1-TOTAL: init:  read_potrho");
    if(ct.runflag <113)
    {
        read_potrho_LCR (vh, vxc, rho);
        if (pct.imgpe == 0) printf ("completed: read_potrho_LCR \n");
    }
    else
    {    
        read_rho_and_pot (ct.infile, vh, vxc, vh_old, vxc_old, rho);
        if (pct.imgpe == 0) printf ("completed: read_rho_and_pot \n");
    }
    delete(RT3);

    if(ct.runflag == 300) 
    {
        plane_average_rho(rho); 
//        exit(0);
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
    MPI_Barrier(pct.img_comm);


    for (level = 0; level < ct.eig_parm.levels + 1; level++)
        make_mask_grid_state (level, states);
    // Initialize some commonly used plans
    FftInitPlans();

    RmgTimer *RT4 = new RmgTimer("1-TOTAL: init:  psp");
    /* Initialize the radial potential stuff */
    InitPseudo(ControlMap);

    /* Initialize the radial qfunction stuff */
    InitQfunct(ControlMap);

    /* Initialize symmetry stuff */
    //init_sym ();

    /* Initialize the nuclear local potential and the compensating charges */
//    init_nuc (vnuc, rhoc, rhocore);
    double *dum_array = NULL;
    InitLocalObject (vnuc, dum_array, ATOMIC_LOCAL_PP, false);
    InitLocalObject (rhoc, dum_array, ATOMIC_RHOCOMP, false);
    InitLocalObject (rhocore, dum_array, ATOMIC_RHOCORE, false);



    delete(RT4);

    init_ext (vext, ct.gbias_begin, ct.gbias_end, ct.BT, ct.gate_bias);

    write_rho_x (rho, "rho_init1");

    RmgTimer *RT5 = new RmgTimer("1-TOTAL: init:  non-local");
    /* Initialize Non-local operators */
    init_nl_xyz ();
    get_ion_orbit_overlap_nl (states);

    GetNlop_on ();
    init_nonlocal_comm ();
    InitNonlocalComm ();

    /* Initialize qfuction in Cartesin coordinates */
    GetQI ();

    /* Get the qqq */
    get_qqq ();
    delete(RT5);



    if (pct.imgpe == 0) printf ("completed: initnegf \n");

}                               /* end init */

