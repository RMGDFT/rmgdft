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
#include "LocalObject.h"
#include "Kbpsi.h"

void InitON(double * vh, double * rho, double *rho_oppo,  double * rhocore, double * rhoc,
          STATE * states, STATE * states1, double * vnuc, double * vxc, double * vh_old, 
          double * vxc_old, std::unordered_map<std::string, InputKey *>& ControlMap)
{

    int ic, idx, ion;
    int level;
    int ione = 1;
    char newname[MAX_PATH + 20];


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
            Atoms[ion].icrds[ic] = Atoms[ion].crds[ic];
            Atoms[ion].ixtal[ic] = Atoms[ion].xtal[ic];
        }

    if (pct.imgpe == 0)
    {
        write_header();
    }

    MPI_Barrier(pct.img_comm);


    RmgTimer *RT = new RmgTimer("1-TOTAL: init: state_init");

    /* allocate memory for wave functions states.psiR and psiI */
    init_state_size(states);

    state_corner_xyz(states);




    int size = (ct.state_end - ct.state_begin) * ct.num_states;

    state_overlap_or_not = new char[size];
    //my_malloc( state_overlap_or_not, size,  char);

    is_state_overlap(states, state_overlap_or_not);


    get_orbit_overlap_region(states);
    GetOrbitalPairs(states);

    delete(RT);

    RmgTimer *RT1 = new RmgTimer("1-TOTAL: init: init_comm");

    init_comm(states);

    init_comm_res(states);
    delete(RT1);
    RmgTimer *RT2 = new RmgTimer("1-TOTAL: init: init_nuc");
    AllocatePsi(states, states1);

    duplicate_states_info(states, states1);
    duplicate_states_info(states, states_tem);

    MPI_Barrier(pct.img_comm);

    {
        int *ixmin, *iymin, *izmin, *dimx, *dimy, *dimz;
        ixmin = new int[3*ct.num_states];
        iymin = ixmin + ct.num_states;
        izmin = iymin + ct.num_states;
        dimx = new int[3*ct.num_states];
        dimy = dimx + ct.num_states;
        dimz = dimy + ct.num_states;

        for(int st = 0; st < ct.num_states; st++)
        {
            ixmin[st] = states[st].ixmin;
            iymin[st] = states[st].iymin;
            izmin[st] = states[st].izmin;
            dimx[st] = states[st].orbit_nx;
            dimy[st] = states[st].orbit_ny;
            dimz[st] = states[st].orbit_nz;

        }


        int density = 1;
        LocalOrbital = new LocalObject<double>(ct.num_states, ixmin, iymin, izmin,
                dimx, dimy, dimz, 0, Rmg_G, density, pct.grid_comm);
        H_LocalOrbital = new LocalObject<double>(ct.num_states, ixmin, iymin, izmin,
                dimx, dimy, dimz, 0, Rmg_G, density, pct.grid_comm);
        delete [] ixmin;
        delete [] dimx;
//        LocalOrbital->ReadOrbitals(std::string(ct.infile), Rmg_G);
    }

    allocate_masks(states);

    for (level = 0; level < ct.eig_parm.levels + 1; level++)
        make_mask_grid_state(level, states);

    // Initialize some commonly used plans
    FftInitPlans();
    /* Initialize the radial potential stuff */
    InitPseudo(ControlMap);

    /* Initialize symmetry stuff */
    //  init_sym();

    /* Initialize the nuclear local potential and the compensating charges */
    //    init_nuc(vnuc, rhoc, rhocore);

    pct.loc_ions_list = new int[ct.num_ions];
    double *dum_array = NULL;
    InitLocalObject (vnuc, dum_array, ATOMIC_LOCAL_PP, false);
    InitLocalObject (rhoc, dum_array, ATOMIC_RHOCOMP, false);
    InitLocalObject (rhocore, dum_array, ATOMIC_RHOCORE, false);


    /* Initialize Non-local operators */
    init_nl_xyz();
    get_ion_orbit_overlap_nl(states);

    GetNlop_on();

    {
        int *ixmin, *iymin, *izmin, *dimx, *dimy, *dimz;
        int tot_prj = 0;
        for (int ion=0; ion < ct.num_ions; ion++)
        {
            tot_prj += Species[Atoms[ion].species].num_projectors;
        }
        ixmin = new int[3*tot_prj];
        iymin = ixmin + tot_prj;
        izmin = iymin + tot_prj;
        dimx = new int[3*tot_prj];
        dimy = dimx + tot_prj;
        dimz = dimy + tot_prj;
        int proj_count = 0;
        int *proj_per_ion = new int[ct.num_ions];
        for (int ion=0; ion < ct.num_ions; ion++)
        {
            ION *iptr = &Atoms[ion];
            SPECIES *sp = &Species[iptr->species];
            proj_per_ion[ion] = sp->num_projectors;
            for (int ip = 0; ip < sp->num_projectors; ip++)
            {
               ixmin[proj_count] = iptr->ixstart;
               iymin[proj_count] = iptr->iystart;
               izmin[proj_count] = iptr->izstart;
                dimx[proj_count] = sp->nldim;
                dimy[proj_count] = sp->nldim;
                dimz[proj_count] = sp->nldim;
                proj_count++;
            }
        }

        int density = 1;
        LocalProj = new LocalObject<double>(tot_prj, ixmin, iymin, izmin,
                dimx, dimy, dimz, 0, Rmg_G, density, pct.grid_comm);
        delete [] ixmin;
        delete [] dimx;
        LocalProj->ReadProjectors(ct.num_ions, ct.max_nlpoints, proj_per_ion, Rmg_G);
        delete [] proj_per_ion;
        Kbpsi_mat = new double[LocalProj->num_tot * LocalOrbital->num_tot]; 
    }

    //if(ct.num_ldaU_ions > 0)
    {
        int tot_orbitals_ldaU = 0;
        for (int ion=0; ion < ct.num_ions; ion++)
        {
            ION *iptr = &Atoms[ion];
            SPECIES *sp = &Species[iptr->species];
            for (int ip = 0; ip < sp->num_orbitals; ip++)
            {
                // This ranges over all orbitals including the m-dependence
                if(sp->awave_is_ldaU[ip])
                {

                    tot_orbitals_ldaU++;
                }
            }
        }

        int *ixmin, *iymin, *izmin, *dimx, *dimy, *dimz;
        int density = 1;
        LocalAtomicOrbital = new LocalObject<double>(tot_orbitals_ldaU, ixmin, iymin, izmin,
                dimx, dimy, dimz, 1, Rmg_G, density, pct.grid_comm);
        LocalAtomicOrbital->GetAtomicOrbitals(ct.num_ions, Rmg_G);

    }


    MPI_Barrier(pct.img_comm);
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
            if(ct.spin_flag)
                sprintf (newname, "%s_spin%d", ct.infile, pct.spinpe);
            else
                sprintf (newname, "%s", ct.infile);

            read_data(newname, vh, vxc, vh_old, vxc_old, rho, vh_corr, states);
            pack_vhstod(vh, ct.vh_ext, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), ct.boundaryflag);
            break;

    }

    if(ct.ON_read_from_RMG)
        ReadDataFromRMG(ct.infile_ON_from_RMG, vh, rho, vxc);

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
            RT1 = new RmgTimer("2-Init: exchange/correlation");
            Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
            F->v_xc(rho, rhocore, ct.XC, ct.vtxc, vxc, ct.spin_flag );
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


    MPI_Barrier(pct.img_comm);
    fflush(NULL);


#if	DEBUG
    print_state_sum(states);
    print_status(states, vh, vxc, vnuc, vh_old, "before leaving init.c  ");
    print_state(&states[0]);
    print_sum(pct.psi_size, states[ct.state_begin].psiR, "init.c states sum ");
#endif

    /* some utilities, used in debuging */


}

