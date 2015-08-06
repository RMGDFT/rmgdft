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
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"
#define min(a,b) (((a)>(b)) ? (b) : (a))

void init_wf_atom(STATE *);
void init_wf_lcao(STATE *);

void init(double * vh, double * rho, double *rho_oppo,  double * rhocore, double * rhoc,
          STATE * states, STATE * states1, double * vnuc, double * vxc, double * vh_old, double * vxc_old)
{

    int ic, idx, ion;
    int level;

   int  gridpe = pct.gridpe;

    /* initialize the lattice basis vectors */

    ct.psi_nbasis = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();

    ct.psi_fnbasis = get_FNX_GRID() * get_FNY_GRID() * get_FNZ_GRID();


    my_malloc_init( ct.energies, ct.max_scf_steps, double );

    ct.states = states;

    init_parameter(states);


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

    /* Initialize the mehrstellen weights */
    get_mehr();


     void *RT = BeginRmgTimer("1-TOTAL: init: state_init");

    /* allocate memory for wave functions states.psiR and psiI */
    init_state_size(states);

    state_corner_xyz(states);




    int size = (ct.state_end - ct.state_begin) * ct.num_states;

    my_malloc( state_overlap_or_not, size,  char);


    is_state_overlap(states, state_overlap_or_not);


    get_orbit_overlap_region(states);
    
    EndRmgTimer(RT);

    void *RT1 = BeginRmgTimer("1-TOTAL: init: init_comm");

    init_comm(states);

    init_comm_res(states);
    EndRmgTimer(RT1);

    void *RT2 = BeginRmgTimer("1-TOTAL: init: init_nuc");
    allocate_psi(states, states1);

    duplicate_states_info(states, states1);
    duplicate_states_info(states, states_tem);

    my_barrier();


    allocate_masks(states);

    for (level = 0; level < ct.eig_parm.levels + 1; level++)
        make_mask_grid_state(level, states);

    /* Initialize the radial potential stuff */
    init_kbr();

    /* Initialize symmetry stuff */
   //  init_sym();

    /* Initialize the nuclear local potential and the compensating charges */
    init_nuc(vnuc, rhoc, rhocore);



    /* Initialize Non-local operators */
    init_nl_xyz();
    get_ion_orbit_overlap_nl(states);

    get_nlop();

    my_barrier();
    EndRmgTimer(RT2);

    void *RT3 = BeginRmgTimer("1-TOTAL: init: init_commi_nonlo");

    init_nonlocal_comm();
    EndRmgTimer(RT3);

    fflush(NULL);

    /* Initialize qfuction in Cartesin coordinates */
    void *RT4 = BeginRmgTimer("1-TOTAL: init: init_qfunc");
    init_qfunct();
    EndRmgTimer(RT4);
    void *RT5 = BeginRmgTimer("1-TOTAL: init: init_QI");

    fflush(NULL);
    get_QI();
    EndRmgTimer(RT5);
    void *RT6 = BeginRmgTimer("1-TOTAL: init: init_qqq");

    /* Get the qqq */
    get_qqq();
    EndRmgTimer(RT6);
    fflush(NULL);

    for (idx = 0; idx < get_FP0_BASIS(); idx++) vh[idx] = ZERO;
    double fac;
    fac = 1.0 - ct.spin_flag * 2.0/3.0 /(1.0 + pct.spinpe);
    if(ct.init_equal_density_flag) fac = 0.5; 
    switch(ct.runflag)
    {
        case 0:
            init_wf(states);

            for (idx = 0; idx < get_FP0_BASIS(); idx++)
            {
                rho[idx] = rhoc[idx] * fac;
                rho_oppo[idx] = rhoc[idx] * (1.0 - fac);
        
            }
            break;
        case LCAO_START:
            init_wf_lcao(states);
            lcao_get_rho(rho);

            break;
        case INIT_FIREBALL:
            init_wf_atom(states);
            init_rho_atom(rho);

            for (idx = 0; idx < get_FP0_BASIS(); idx++)
            {
                rho_oppo[idx] = rho[idx] * (1.0 - fac);
                rho[idx] = rho[idx] * fac;
            }
            break;

        case INIT_GAUSSIAN:
            init_wf_gaussian(states);
            for (idx = 0; idx < get_FP0_BASIS(); idx++)
            {
                rho[idx] = rhoc[idx] * fac;
                rho_oppo[idx] = rhoc[idx] * (1.0 - fac);
            }

            break;

        case 1:
        case 5:
            read_data(ct.infile, vh, vxc, vh_old, vxc_old, rho, states);
            pack_vhstod(vh, ct.vh_ext, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), ct.boundaryflag);
            break;

    }

    switch(ct.runflag)
    {
        case 0:
        case LCAO_START:
        case INIT_FIREBALL:
        case INIT_GAUSSIAN:
            get_vxc(rho, rho_oppo, rhocore, vxc);
            pack_vhstod(vh, ct.vh_ext, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), ct.boundaryflag);

            for (idx = 0; idx < get_FP0_BASIS(); idx++)
                rho_tot[idx] = rho[idx] + rho_oppo[idx];

            get_vh (rho_tot, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, 0.0, ct.boundaryflag);
            for (idx = 0; idx < get_FP0_BASIS(); idx++)
                vh_old[idx] = vh[idx];
            for (idx = 0; idx < get_FP0_BASIS(); idx++)
                vxc_old[idx] = vxc[idx];
    }


    void *RT8 = BeginRmgTimer("1-TOTAL: init: init_ddd");
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx];

    get_vtot_psi(vtot_c, vtot, get_FG_RATIO());
    get_ddd(vtot);
    EndRmgTimer(RT8);


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
