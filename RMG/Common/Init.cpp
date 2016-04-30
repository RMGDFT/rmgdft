/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "Functional.h"
#include "GpuAlloc.h"
#include "../Headers/prototypes.h"
#include "ErrorFuncs.h"
#include "RmgException.h"
#include "Functional.h"
#include "Solvers.h"
#include "pfft.h"
#include "RmgParallelFft.h"

static void init_alloc_nonloc_mem (void);


template <typename OrbitalType> void Init (double * vh, double * rho, double * rho_oppo, double * rhocore, double * rhoc,
           double * vnuc, double * vxc,  Kpoint<OrbitalType> **Kptr);

// Instantiate gamma and non-gamma versions
template void Init<double>(double*, double*, double*, double*, double*, double*, double*, Kpoint<double>**);
template void Init<std::complex<double> >(double*, double*, double*, double*, double*, double*, double*, Kpoint<std::complex <double> >**);

template <typename OrbitalType> void Init (double * vh, double * rho, double * rho_oppo, double * rhocore, double * rhoc,
           double * vnuc, double * vxc,  Kpoint<OrbitalType> **Kptr)
{

    RmgTimer RT0("Init");
    int kpt, kpt1, ic, idx, state, ion, st1, P0_BASIS, FP0_BASIS;
    int species;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID;

    SPECIES *sp;
    OrbitalType *rptr = NULL, *nv, *ns, *Bns = NULL;
    double *vtot, *rho_tot;
    double time2, fac;

    nv = (OrbitalType *)pct.nv;

    P0_BASIS =  Rmg_G->get_P0_BASIS(1);
    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    FPX0_GRID = Rmg_G->get_PX0_GRID(Rmg_G->get_default_FG_RATIO());
    FPY0_GRID = Rmg_G->get_PY0_GRID(Rmg_G->get_default_FG_RATIO());
    FPZ0_GRID = Rmg_G->get_PZ0_GRID(Rmg_G->get_default_FG_RATIO());

    ct.fftw_wisdom_setup = 0;

    ct.total_scf_steps = 0;
    ct.md_steps = 0;

    ct.psi_nbasis = Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1);
    ct.psi_fnbasis = Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO) * Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO) * Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO);

#if SCALAPACK_LIBS
    if(ct.subdiag_driver == SUBDIAG_SCALAPACK) {

        /*Initialize ScaLapack, get context */
        sl_init (&pct.ictxt, ct.num_states);

        /*Set desca, variable used in ScaLapack functions */
        set_desca (pct.desca, &pct.ictxt, ct.num_states);

    }
#endif

    // Initialize some commonly used plans for our parallel ffts
    FftInitPlans();

    /* Allocate storage for non-local projectors */
    pct.newsintR_local = NULL;
    pct.oldsintR_local = NULL;


    /* Set hartree boundary condition stuff */
    ct.vh_pxgrid = FPX0_GRID;
    ct.vh_pygrid = FPY0_GRID;
    ct.vh_pzgrid = FPZ0_GRID;

    ct.vh_pbasis = ct.vh_pxgrid * ct.vh_pygrid * ct.vh_pzgrid;
    ct.vh_ext = new double[ct.vh_pbasis];
    
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        vh[idx] = 0.0;
        ct.vh_ext[idx] = 0.0;
    }



    /* Get the crystal or cartesian coordinates of the ions */
    init_pos ();


    ct.hmaxgrid = Rmg_L.get_xside() * Rmg_G->get_hxgrid(1);
    if (Rmg_L.get_yside() * Rmg_G->get_hygrid(1) > ct.hmaxgrid)
        ct.hmaxgrid = Rmg_L.get_yside() * Rmg_G->get_hygrid(1);
    if (Rmg_L.get_zside() * Rmg_G->get_hzgrid(1) > ct.hmaxgrid)
        ct.hmaxgrid = Rmg_L.get_zside() * Rmg_G->get_hzgrid(1);

    ct.hmingrid = Rmg_L.get_xside() * Rmg_G->get_hxgrid(1);
    if (Rmg_L.get_yside() * Rmg_G->get_hygrid(1) < ct.hmingrid)
        ct.hmingrid = Rmg_L.get_yside() * Rmg_G->get_hygrid(1);
    if (Rmg_L.get_zside() * Rmg_G->get_hzgrid(1) < ct.hmingrid)
        ct.hmingrid = Rmg_L.get_zside() * Rmg_G->get_hzgrid(1);


    if ((ct.hmaxgrid / ct.hmingrid) > 1.1)
    {
        if (pct.imgpe == 0)
        {
            printf ("hxgrid = %7.5f\n", Rmg_G->get_hxgrid(1));
            printf ("hygrid = %7.5f\n", Rmg_G->get_hygrid(1));
            printf ("hzgrid = %7.5f\n", Rmg_G->get_hzgrid(1));
        }
        rmg_error_handler (__FILE__, __LINE__, "Anisotropy too large");
    }


    /* Set discretization array */
    ct.xcstart = 0.0;
    ct.ycstart = 0.0;
    ct.zcstart = 0.0;

    // Compute some buffer sizes. Have to deal with the interaction of a couple of different options here
    //
    //    If an LCAO start is selected we need to allocate sufficient memory for the initial set of orbitals
    //    which may be larger than the number of orbitals used for the rest of the run.
    //    If potential acceleration is selected we need another buffer of size (run_states * P0_BASIS). If
    //    we allocate this in a single large buffer of 2*run_states*P0_BASIS it will probably be enough for
    //    the initialization even if an LCAO start is selected. Finally if Davidson diagonalization is
    //    requested we need to allocate memory for the expanded basis including the diagonalization arrays.

    bool potential_acceleration = ((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0));
    


    /* Set state pointers and initialize state data */
#if GPU_ENABLED

    cudaError_t custat;

    // Figure out how much memory space to reserve on the GPU
    // 3 blocks of RMG_CUBLASXT_BLOCKSIZE*ct.max_states for diagonalization arrays
    size_t gpu_bufsize, t1;
    t1 = RMG_CUBLASXT_BLOCKSIZE * ct.max_states * sizeof(OrbitalType);
    gpu_bufsize = 3 * t1;
#if MAGMA_LIBS
    gpu_bufsize += t1;
#endif

    InitGpuMalloc(gpu_bufsize);

    // Next is page locked memory for transferring data back and forth
    int n_win = 0;
    if(ct.use_folded_spectrum) {
        double r_width = ct.folded_spectrum_width;
        double t1 = (double)ct.num_states;
        n_win = (int)(r_width * t1) + 1;
    }

    size_t gpu_hostbufsize;
    gpu_hostbufsize = 2 * ct.max_states * ct.max_states * sizeof(OrbitalType) + 
                      3 * ct.max_states * std::max(ct.max_states, P0_BASIS) * sizeof(OrbitalType) +
                      n_win * n_win * sizeof(OrbitalType);

    InitGpuMallocHost(gpu_hostbufsize);

    // Wavefunctions are actually stored here
    custat = cudaMallocHost((void **)&rptr, (ct.num_kpts * ct.max_states * P0_BASIS + 1024) * sizeof(OrbitalType));
    RmgCudaError(__FILE__, __LINE__, custat, "cudaMallocHost failed");
    custat = cudaMallocHost((void **)&nv, ct.non_local_block_size * P0_BASIS * sizeof(OrbitalType));
    RmgCudaError(__FILE__, __LINE__, custat, "cudaMallocHost failed");
    custat = cudaMallocHost((void **)&ns, ct.max_states * P0_BASIS * sizeof(OrbitalType));
    RmgCudaError(__FILE__, __LINE__, custat, "cudaMallocHost failed");
    if(!ct.norm_conserving_pp) {
        custat = cudaMallocHost((void **)&Bns, ct.non_local_block_size * P0_BASIS * sizeof(OrbitalType));
        RmgCudaError(__FILE__, __LINE__, custat, "cudaMallocHost failed");
        pct.Bns = (double *)Bns;
    }
#else
    // Wavefunctions are actually stored here
    rptr = new OrbitalType[ct.num_kpts * ct.max_states * P0_BASIS + 1024];
    nv = new OrbitalType[ct.non_local_block_size * P0_BASIS]();
    ns = new OrbitalType[ct.max_states * P0_BASIS]();
    if(!ct.norm_conserving_pp) {
        Bns = new OrbitalType[ct.non_local_block_size * P0_BASIS]();
        pct.Bns = (double *)Bns;
    }
#endif
    pct.nv = (double *)nv;
    pct.ns = (double *)ns;


    kpt1 = ct.num_kpts;
//    if (Verify ("calculation_mode", "Band Structure Only", Kptr[0]->ControlMap))
//        kpt1 = 1;
    for (kpt = 0; kpt < kpt1; kpt++)
    {

        Kptr[kpt]->set_pool(rptr);
        Kptr[kpt]->nv = nv;
        Kptr[kpt]->ns = ns;

        if(!ct.norm_conserving_pp) Kptr[kpt]->Bns = Bns;

        for (st1 = 0; st1 < ct.max_states; st1++)
        {
            Kptr[kpt]->Kstates[st1].kidx = kpt;
            Kptr[kpt]->Kstates[st1].psi = rptr;
            Kptr[kpt]->Kstates[st1].vxc = vxc;
            Kptr[kpt]->Kstates[st1].vh = vh;
            Kptr[kpt]->Kstates[st1].vnuc = vnuc;
            Kptr[kpt]->Kstates[st1].pbasis = P0_BASIS;
            Kptr[kpt]->Kstates[st1].istate = st1;
            rptr +=P0_BASIS;
        }
    }


    //Dprintf ("If not an initial run read data from files");
    if ((ct.runflag == RESTART) || (ct.forceflag == BAND_STRUCTURE) )
    {
        int THREADS_PER_NODE = ct.THREADS_PER_NODE;
        ReadData (ct.infile, vh, rho, vxc, Kptr);
        ct.THREADS_PER_NODE = THREADS_PER_NODE;

        /*For spin polarized calculation we need to get opposite charge density, eigenvalues and occupancies*/
        if (ct.spin_flag)
        {
            get_rho_oppo (rho, rho_oppo);
            GetOppositeEigvals (Kptr);
            GetOppositeOccupancies (Kptr);
        }
    }
    else 
    {
        /* Set the initial hartree potential to a constant */
        for (idx = 0; idx < FP0_BASIS; idx++)
            vh[idx] = 0.0;
    }

    /* Set initial states to random start */

    if (ct.runflag == LCAO_START || (ct.forceflag == BAND_STRUCTURE))
    {
        RmgTimer *RT2 = new RmgTimer("Init: LcaoGetPsi");
        for(kpt = 0; kpt < ct.num_kpts; kpt++) {
            LcaoGetPsi(Kptr[kpt]->Kstates);
        }
        delete(RT2);
    }

    if (ct.runflag == RANDOM_START)
    {
        for (kpt = 0; kpt < ct.num_kpts; kpt++)
            //init_wf (&Kptr[0]->kstates[kpt * ct.num_states]);
            Kptr[kpt]->random_init();
    }


    /* Set initial ionic coordinates to the current ones. */
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        for (ic = 0; ic < 3; ic++)
        {
            ct.ions[ion].icrds[ic] = ct.ions[ion].crds[ic];
            ct.ions[ion].ixtal[ic] = ct.ions[ion].xtal[ic];
        }
    }




    //Dprintf ("Initialize the radial potential stuff");
    RmgTimer *RT1 = new RmgTimer("Init: radial potentials");
    InitPseudo (Kptr[0]->ControlMap);
    delete(RT1);

    /* Initialize symmetry stuff */
    if(!ct.is_gamma)
        init_sym ();

    //Dprintf ("Allocate memory for arrays related to nonlocal PP");
    init_alloc_nonloc_mem ();

    /*Set max_nldim */
    ct.max_nldim = 0;
    for (species = 0; species < ct.num_species; species++)
    {
        /* Get species type */
        sp = &ct.sp[species];

        if (sp->nldim > ct.max_nldim)
            ct.max_nldim = sp->nldim;
    }

    time2 = my_crtc ();

    rmg_printf ("\n\n init: Starting FFTW initialization ...");

    fflush (NULL);


    /*Do forward transform for each species and store results on the coarse grid */
    RT1 = new RmgTimer("Init: weights");
    init_weight ();
    /*The same for derivative of beta */
    init_derweight ();
    delete(RT1);

    rmg_printf ("\n init: FFTW initialization finished, it took %.1f s", my_crtc () - time2);
    fflush (NULL);


    /* Initialize the qfunction stuff */
    RT1 = new RmgTimer("Init: qfunct");
    InitQfunct(Kptr[0]->ControlMap);
    delete(RT1);

    /* Update items that change when the ionic coordinates change */
    RT1 = new RmgTimer("Init: ionic potentials");
    ReinitIonicPotentials (Kptr, vnuc, rhocore, rhoc);
    delete(RT1);


    // Normalize orbitals if not an initial run
    if (ct.runflag != RESTART) /* Initial run */
    {
        RmgTimer *RT2 = new RmgTimer("Init: normalization");
        for (kpt = 0; kpt < ct.num_kpts; kpt++)
        {
            for (state = 0; state < ct.num_states; state++)
            {
                Kptr[kpt]->Kstates[state].normalize(Kptr[kpt]->Kstates[state].psi, state);
            }
        }
        delete(RT2);
    }

    //mix_betaxpsi(0);

    //Dprintf ("Set the initial density to be equal to the compensating charges");

    /*For random start, use charge density equal to compensating charge */
    if (ct.runflag == RANDOM_START)
    {
        if (ct.spin_flag)
        {   
            fac = (2.0 - ct.init_equal_density_flag) / (3.0 - ct.init_equal_density_flag);
            for (idx = 0; idx < FP0_BASIS; idx++)
            {

                if (pct.spinpe == 0)
                {				    
                    rho[idx] = rhoc[idx] * fac;
                    rho_oppo[idx] = rhoc[idx] * (1.0 - fac);
                }
                else         /* if pct.spinpe = 1  */
                {
                    rho_oppo[idx] = rhoc[idx] * fac;
                    rho[idx] = rhoc[idx] * (1.0 - fac);
                }
            }

        }

        else
        {
            for (idx = 0; idx < FP0_BASIS; idx++)
                rho[idx] = rhoc[idx];
        }
    }

    if (ct.runflag == LCAO_START && (ct.forceflag != BAND_STRUCTURE))
        lcao_get_rho(rho);


    ct.rms = 0.0;

    //Dprintf ("Generate initial vxc potential and hartree potential");
    pack_vhstod (vh, ct.vh_ext, FPX0_GRID, FPY0_GRID, FPZ0_GRID, ct.boundaryflag);



    //Dprintf ("Condition of run flag is %d", ct.runflag);
    /*If not a restart, get vxc and vh, for restart these quantities should read from the restart file*/
    if (ct.runflag != RESTART)
    {
        //get_vxc (rho, rho_oppo, rhocore, vxc);
        double etxc, vtxc;
        Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
        F->v_xc(rho, rhocore, etxc, vtxc, vxc, ct.spin_flag );
        delete F;

        if (ct.spin_flag)
        {
            rho_tot = new double[FP0_BASIS];
            for (idx = 0; idx < FP0_BASIS; idx++)
                rho_tot[idx] = rho[idx] + rho_oppo[idx];

            get_vh (rho_tot, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, 0.0, ct.boundaryflag);
            delete [] rho_tot;
        }
        else
            get_vh (rho, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, 0.0, ct.boundaryflag);

    }

    // Generate initial Betaxpsi
    for(kpt = 0;kpt < ct.num_kpts;kpt++) {
        Betaxpsi (Kptr[kpt], 0, Kptr[kpt]->nstates, Kptr[kpt]->newsint_local, Kptr[kpt]->nl_weight);
        Kptr[kpt]->mix_betaxpsi(0);
    }

    // If not a restart and diagonalization is requested do a subspace diagonalization otherwise orthogonalize
    if(ct.runflag != RESTART) {

        if (ct.initdiag)
        {
            /*dnmI has to be stup before calling subdiag */
            vtot = new double[FP0_BASIS];
            double *vtot_psi = new double[P0_BASIS];

            for (idx = 0; idx < FP0_BASIS; idx++)
                vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];

            /*Generate the Dnm_I */
            get_ddd (vtot);

            // Transfer vtot from the fine grid to the wavefunction grid for Subdiag
            GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);

            /*Now we can do subspace diagonalization */
            double *new_rho=new double[FP0_BASIS];
            for(kpt = 0;kpt < ct.num_kpts;kpt++) {
                Subdiag (Kptr[kpt], vtot_psi, ct.subdiag_driver);
                Betaxpsi (Kptr[kpt], 0, Kptr[kpt]->nstates, Kptr[kpt]->newsint_local, Kptr[kpt]->nl_weight);
                Kptr[kpt]->mix_betaxpsi(0);
            }
            // Get new density 
            GetNewRho(Kptr, new_rho);
            MixRho(new_rho, rho, rhocore, vh, vh, rhoc, Kptr[0]->ControlMap);
            delete [] new_rho;

            /*Release vtot memory */
            delete [] vtot_psi;
            delete [] vtot;


        }                           /*end if(ct.initdiag) */
        else
        {

            for (kpt = 0; kpt < ct.num_kpts; kpt++) {
                Kptr[kpt]->orthogonalize(Kptr[kpt]->orbital_storage);
            }

        }

    }


    ct.num_states = ct.run_states;
    for (kpt = 0; kpt < ct.num_kpts; kpt++) {
        Kptr[kpt]->nstates = ct.run_states;
        Kptr[kpt]->dvh_skip = 8;
        // Set up potential acceleration arrays if required
        if(potential_acceleration) {
            if(ct.run_states <= 256) Kptr[kpt]->dvh_skip = 4;
            if(ct.run_states <= 128) Kptr[kpt]->dvh_skip = 2;
            if(ct.run_states <= 64) Kptr[kpt]->dvh_skip = 1;
            Kptr[kpt]->ndvh = ct.run_states / Kptr[kpt]->dvh_skip + 1;
            Kptr[kpt]->dvh = new double[ Kptr[kpt]->ndvh * P0_BASIS ];
        }
    }

    if (Verify ("start_mode","LCAO Start", Kptr[0]->ControlMap) && (ct.num_states != ct.init_states)) {

#if SCALAPACK_LIBS
    if(ct.subdiag_driver == SUBDIAG_SCALAPACK) {
        // In case matrix sizes changed
//        if (pct.scalapack_pe) {
//            sl_exit(pct.ictxt, 0);
//        }
//        sl_init (&pct.ictxt, ct.num_states);
//        set_desca (pct.desca, &pct.ictxt, ct.num_states);
    }
#endif

    }

}                               /* end init */




static void init_alloc_nonloc_mem (void)
{
    int ion;


    pct.nlindex = new int *[ct.num_ions];
    pct.Qindex = new int *[ct.num_ions];
    pct.idxflag = new int *[ct.num_ions];
    pct.Qdvec = new int *[ct.num_ions];

    pct.idxptrlen = new int [ct.num_ions];
    pct.Qidxptrlen = new int [ct.num_ions];
    pct.lptrlen = new int [ct.num_ions];

    pct.phaseptr = new double *[ct.num_ions];
    pct.augfunc = new double *[ct.num_ions];
    pct.dnmI = new double *[ct.num_ions];
    pct.qqq = new double *[ct.num_ions];


    /*Initialize pointer arrays to NULL */
    pct.weight = NULL;
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        ION *iptr = &ct.ions[ion];

        iptr->fftw_phase_sin = NULL;
        iptr->fftw_phase_cos = NULL;

        pct.idxptrlen[ion] = 0;
        pct.Qidxptrlen[ion] = 0;
        pct.lptrlen[ion] = 0;

        pct.nlindex[ion] = NULL;
        pct.Qindex[ion] = NULL;
        pct.idxflag[ion] = NULL;
        pct.Qdvec[ion] = NULL;

        pct.phaseptr[ion] = NULL;
        pct.augfunc[ion] = NULL;
        pct.dnmI[ion] = NULL;
        pct.qqq[ion] = NULL;

    }                           /*end for(ion=0; ion<ct.num_ions; ion++) */

}                               /*end init_alloc_nonloc_mem */

/******/
