#include "negf_prototypes.h"
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
#include "LocalObject.h"
#include "GpuAlloc.h"




void is_state_overlap (STATE *, char *);

void InitNegf (double * vh, double * rho, double * rhocore, double * rhoc, double * rho_tf,
                STATE * states, STATE * states1, double * vnuc, double * vext, double * vxc, double * vh_old,
                double * vxc_old, std::unordered_map<std::string, InputKey *>& ControlMap)
{

    int ic, idx, ion;
    double tem;
    int st1, iprobe, i;

    ct.psi_nbasis = Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1);

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
            Atoms[ion].icrds[ic] = Atoms[ion].crds[ic];
            Atoms[ion].ixtal[ic] = Atoms[ion].xtal[ic];
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

    init_state_size (states);

    state_corner_xyz (states);

    MPI_Barrier(pct.img_comm);

    pmo_init();


    int *ixmin, *iymin, *izmin, *dimx, *dimy, *dimz;
    ixmin = new int[6*ct.num_states];
    iymin = ixmin + ct.num_states;
    izmin = ixmin + 2*ct.num_states;
    dimx  = ixmin + 3*ct.num_states;
    dimy  = ixmin + 4*ct.num_states;
    dimz  = ixmin + 5*ct.num_states;

    int density = 1;
    for(int st = 0; st < ct.num_states; st++)
    {

        ixmin[st] = states[st].ixmin;
        iymin[st] = states[st].iymin;
        izmin[st] = states[st].izmin;
        dimx[st] = states[st].orbit_nx;
        dimy[st] = states[st].orbit_ny;
        dimz[st] = states[st].orbit_nz;

    }

    LocalOrbital = new LocalObject<double>(ct.num_states, ixmin, iymin, izmin,
            dimx, dimy, dimz, 0, *Rmg_G, density, pct.grid_comm);
    H_LocalOrbital = new LocalObject<double>(ct.num_states, ixmin, iymin, izmin,
            dimx, dimy, dimz, 0, *Rmg_G, density, pct.grid_comm);
    LocalOrbital->ReAssign(*Rmg_G);
    H_LocalOrbital->ReAssign(*Rmg_G);


    delete [] ixmin;


    RmgTimer *RT1 = new RmgTimer("1-TOTAL: init:  read_orbital");
    ReadInterpolateOrbitals();
    delete(RT1);

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


    // Initialize some commonly used plans
    FftInitPlans();

    RmgTimer *RT4 = new RmgTimer("1-TOTAL: init:  psp");
    /* Initialize the radial potential stuff */
    InitPseudo();

    /* Initialize the radial qfunction stuff */
    InitQfunct();

    /* Initialize symmetry stuff */
    //init_sym ();

    /* Initialize the nuclear local potential and the compensating charges */
    //    init_nuc (vnuc, rhoc, rhocore);
    pct.loc_ions_list = new int[ct.num_ions];
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

    GetNlop_on ();

    /* Initialize qfuction in Cartesin coordinates */
    GetQI ();

    /* Get the qqq */
    get_qqq ();

    int tot_prj = 0;
    for (int ion=0; ion < ct.num_ions; ion++)
    {
        tot_prj += Species[Atoms[ion].species].num_projectors;
    }
    ixmin = new int[6*tot_prj];
    iymin = ixmin + tot_prj;
    izmin = ixmin + 2*tot_prj;
    dimx = ixmin + 3*tot_prj;
    dimy = ixmin + 4*tot_prj;
    dimz = ixmin + 5*tot_prj;
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

    LocalProj = new LocalObject<double>(tot_prj, ixmin, iymin, izmin,
            dimx, dimy, dimz, 0, *Rmg_G, density, pct.grid_comm);
    LocalProj->ReadProjectors(ct.num_ions, ct.max_nlpoints, proj_per_ion, *Rmg_G);

    delete [] ixmin;
    delete [] proj_per_ion;

    size_t size = LocalProj->num_thispe * LocalOrbital->num_thispe * sizeof(double);
    Kbpsi_mat_local = (double *) GpuMallocManaged(size);

    size = LocalProj->num_tot * LocalOrbital->num_tot * sizeof(double);
    Kbpsi_mat = (double *) GpuMallocManaged(size);
    for(int ib = 0; ib < ct.num_blocks; ib++)
        Kbpsi_mat_blocks.push_back(&Kbpsi_mat[ pmo.orb_index[ib] * LocalProj->num_tot]); 

    delete(RT5);

    if (pct.imgpe == 0) printf ("completed: initnegf \n");

}                               /* end init */

