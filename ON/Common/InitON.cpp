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
 *   init_sys.c get_nlop.c init_nuc.c init_wflcao.c get_vxc.c
 *   pack_vhstod.c 
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


#include "blas.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "RmgParallelFft.h"

#include "prototypes_on.h"
#include "init_var.h"
#include "transition.h"
#include "Atomic.h"
#include "Functional.h"

void InitON(double * vh, double * rho, double *rho_oppo,  double * rhocore, double * rhoc,
          STATE * states, STATE * states1, double * vnuc, double * vxc, double * vh_old, 
          double * vxc_old, std::unordered_map<std::string, InputKey *>& ControlMap)
{

    int ic, idx, ion;
    int level;
    int ione = 1;


    /* initialize the lattice basis vectors */

    ct.psi_nbasis = Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1);

    ct.psi_fnbasis = Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO) *
        Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO) *
        Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO) ;

    ct.energies = new double[ct.max_scf_steps];

    ct.states = states;


    init_parameter(states);

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

    if (pct.imgpe == 0)
    {
        write_header();
    }

    my_barrier();


    RmgTimer *RT = new RmgTimer("1-TOTAL: init: state_init");

    /* allocate memory for wave functions states.psiR and psiI */
    init_state_size(states);

    state_corner_xyz(states);




    int size = (ct.state_end - ct.state_begin) * ct.num_states;

    state_overlap_or_not = new char[size];
    //my_malloc( state_overlap_or_not, size,  char);


    is_state_overlap(states, state_overlap_or_not);


    get_orbit_overlap_region(states);

    delete(RT);

    RmgTimer *RT1 = new RmgTimer("1-TOTAL: init: init_comm");

    init_comm(states);

    init_comm_res(states);
    delete(RT1);
    RmgTimer *RT2 = new RmgTimer("1-TOTAL: init: init_nuc");
    allocate_psi(states, states1);

    duplicate_states_info(states, states1);
    duplicate_states_info(states, states_tem);

    my_barrier();


    allocate_masks(states);

    for (level = 0; level < ct.eig_parm.levels + 1; level++)
        make_mask_grid_state(level, states);

    /* Initialize the radial potential stuff */
    InitPseudo(ControlMap);

    /* Initialize symmetry stuff */
    //  init_sym();

    /* Initialize the nuclear local potential and the compensating charges */
//    init_nuc(vnuc, rhoc, rhocore);

    double *dum_array = NULL;
    InitLocalObject (vnuc, dum_array, ATOMIC_LOCAL_PP, false);
    InitLocalObject (rhoc, dum_array, ATOMIC_RHOCOMP, false);
    InitLocalObject (rhocore, dum_array, ATOMIC_RHOCORE, false);


    /* Initialize Non-local operators */
    init_nl_xyz();
    get_ion_orbit_overlap_nl(states);

    GetNlop_on();

    my_barrier();
    delete(RT2);

    RmgTimer *RT3 = new RmgTimer("1-TOTAL: init: init_commi_nonlo");

    init_nonlocal_comm();
    InitNonlocalComm();
    delete(RT3);

    fflush(NULL);

    /* Initialize qfuction in Cartesin coordinates */
    RmgTimer *RT4 = new RmgTimer("1-TOTAL: init: init_qfunc");
    InitQfunct(ControlMap);
    delete(RT4);
    RmgTimer *RT5 = new RmgTimer("1-TOTAL: init: init_QI");

    fflush(NULL);
    GetQI();
    delete(RT5);
    RmgTimer *RT6 = new RmgTimer("1-TOTAL: init: init_qqq");

    /* Get the qqq */
    get_qqq();
    delete(RT6);
    fflush(NULL);

    int np[3];
    ptrdiff_t densgrid[3];
    np[0] = Rmg_G->get_PE_X();
    np[1] = Rmg_G->get_PE_Y();
    np[2] = Rmg_G->get_PE_Z();
    densgrid[0] =  Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO);
    densgrid[1] =  Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO);
    densgrid[2] =  Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO);

    // Initialize some commonly used plans
    FftInitPlans();


    if(ct.dipole_corr[0] + ct.dipole_corr[1] + ct.dipole_corr[2] > 0)
    {
        VhcorrDipoleInit(vh_x, vh_y, vh_z, rhoc);
    }


    int FP0_BASIS = get_FP0_BASIS();
    for (idx = 0; idx < get_FP0_BASIS(); idx++) vh[idx] = ZERO;
    switch(ct.runflag)
    {
        case 0:
            init_wf(states);

            dcopy(&FP0_BASIS, rhoc, &ione, rho, &ione);

            break;
        case LCAO_START:
            init_wf_lcao(states);
            LcaoGetAtomicRho(rho);

            break;
        case INIT_FIREBALL:
            init_wf_atom(states);
            init_rho_atom(rho);

            break;

        case INIT_GAUSSIAN:
            init_wf_gaussian(states);
            dcopy(&FP0_BASIS, rhoc, &ione, rho, &ione);

            break;

        case 1:
        case 5:
        case 6:
            read_data(ct.infile, vh, vxc, vh_old, vxc_old, rho, vh_corr, states);
            pack_vhstod(vh, ct.vh_ext, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), ct.boundaryflag);
            break;

    }

    if(ct.spin_flag) 
        get_rho_oppo (rho,  rho_oppo);
    else
        for (idx = 0; idx < get_FP0_BASIS(); idx++)
        {
            rho_oppo[idx] = 0.0;
        }


    switch(ct.runflag)
    {
        case 0:
        case LCAO_START:
        case INIT_FIREBALL:
        case INIT_GAUSSIAN:
            double tcharge = 0.0;
            for (idx = 0; idx < get_FP0_BASIS(); idx++)
                tcharge += rho[idx];


            ct.tcharge = real_sum_all(tcharge, pct.grid_comm);
            ct.tcharge = real_sum_all(ct.tcharge, pct.spin_comm);


            ct.tcharge *= get_vel_f();

            double t2 = ct.nel / ct.tcharge;
            int iii = get_FP0_BASIS();
            dscal(&iii, &t2, &rho[0], &ione);

            //    get_vxc(rho, rho_oppo, rhocore, vxc);
            double vtxc, etxc;
            RT1 = new RmgTimer("2-Init: exchange/correlation");
            Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
            F->v_xc(rho, rhocore, etxc, vtxc, vxc, ct.spin_flag );
            delete F;
            delete RT1;
            pack_vhstod(vh, ct.vh_ext, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), ct.boundaryflag);

            double rms_target = ct.rms/ct.hartree_rms_ratio;
            VhDriver(rho, rhoc, vh, ct.vh_ext, rms_target);

            for (idx = 0; idx < get_FP0_BASIS(); idx++)
            {
                vh_old[idx] = vh[idx];
                vxc_old[idx] = vxc[idx];
                vh_corr[idx] = 0.0;
            }
    }


    RmgTimer *RT8 = new RmgTimer("1-TOTAL: init: init_ddd");
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx] + vh_corr[idx];

    get_vtot_psi(vtot_c, vtot, get_FG_RATIO());
    get_ddd(vtot);
    delete(RT8);


    my_barrier();
    fflush(NULL);


#if	DEBUG
    print_state_sum(states);
    print_status(states, vh, vxc, vnuc, vh_old, "before leaving init.c  ");
    print_state(&states[0]);
    print_sum(pct.psi_size, states[ct.state_begin].psiR, "init.c states sum ");
#endif

    /* some utilities, used in debuging */


}


