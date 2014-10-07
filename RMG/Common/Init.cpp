/************************** SVN Revision Information **************************
 **    $Id: init.c 2303 2014-05-03 13:14:47Z ebriggs $    **
******************************************************************************/

/****f* QMD-MGDFT/init.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void init(double *vh, double *rho, double *rhocore, double *rhoc, 
 *             STATE *states, double *vnuc, double *vxc)
 *   Basic initialization stuff.
 * INPUTS
 *   rhocore: core charge density for non-linear core corection
 *            it was read in main.c by call read_pseudo.c
 * OUTPUT
 *   vh:  initial Hartree potential 
 *   rho: initial charge density
 *   rhoc: compensating charge density
 *   states: wave functions are initialize randomly or read from data
 *   vnuc:  pseudopotential
 *   vxc: exchange correlation potential calculated from initial rho
 * PARENTS
 *   main.c
 * CHILDREN
 *   Latgen.c recips.c init_pos.c read_data.c init_wf.c init_kbr.c
 *   init_sys.c get_nlop.c init_nuc.c get_vxc.c
 *   pack_vhstod.c get_vh.c
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
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
#include "GpuAlloc.h"
#include "../Headers/prototypes.h"
#include "ErrorFuncs.h"

static void init_alloc_nonloc_mem (void);

extern "C" bool     verify( char *tagname, const void *optvalue );

template <typename OrbitalType> void Init (double * vh, double * rho, double * rho_oppo, double * rhocore, double * rhoc,
           double * vnuc, double * vxc,  Kpoint<OrbitalType> **Kptr);

// Instantiate gamma and non-gamma versions
template void Init<double>(double*, double*, double*, double*, double*, double*, double*, Kpoint<double>**);
template void Init<std::complex<double> >(double*, double*, double*, double*, double*, double*, double*, Kpoint<std::complex <double> >**);

template <typename OrbitalType> void Init (double * vh, double * rho, double * rho_oppo, double * rhocore, double * rhoc,
           double * vnuc, double * vxc,  Kpoint<OrbitalType> **Kptr)
{

    int kpt, kpt1, ic, idx, state, ion, st1, P0_BASIS, FP0_BASIS;
    int species;
    int PX0_GRID, PY0_GRID, PZ0_GRID;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID;

    SPECIES *sp;
    OrbitalType *rptr = NULL, *nv, *ns, *Bns;
    double *vtot, *rho_tot, *rptr1=NULL;
    ION *iptr;
    double time2, fac;

#if GPU_ENABLED
    init_gpu();
#endif

    nv = (OrbitalType *)pct.nv;

    P0_BASIS =  Rmg_G->get_P0_BASIS(1);
    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    PX0_GRID = Rmg_G->get_PX0_GRID(1);
    PY0_GRID = Rmg_G->get_PY0_GRID(1);
    PZ0_GRID = Rmg_G->get_PZ0_GRID(1);
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


    /* Some multigrid parameters */
    // ct.poi_parm.sb_step = 1.0;
    ct.eig_parm.sb_step = 1.0;

    /* Set state pointers and initialize state data */
#if GPU_ENABLED

    cudaError_t custat;

    // Figure out how much memory space to reserve on the GPU
    // 3 blocks of num_states * num_states for diagonalization arrays
    size_t gpu_bufsize, t1;
    t1 = ct.num_states * ct.num_states * sizeof(OrbitalType);
    gpu_bufsize = 4 * t1;
#if MAGMA_LIBS
    gpu_bufsize += t1;
#endif

    // Two buffers for rotating the orbitals. Make sure they are big enough to use
    // as additional diagonalization arrays as well
    t1 = ct.num_states * std::max(ct.num_states, P0_BASIS) * sizeof(OrbitalType);
    gpu_bufsize += 2 * t1;
    // and multiply by 2 just for kicks
    //gpu_bufsize *= 2;
    InitGpuMalloc(gpu_bufsize);

    // Next is page locked memory for transferring data back and forth
    size_t gpu_hostbufsize;
    gpu_hostbufsize = 7 * ct.num_states * ct.num_states * sizeof(OrbitalType);
    InitGpuMallocHost(gpu_hostbufsize);

    // Wavefunctions are actually stored here
    custat = cudaMallocHost((void **)&rptr, (ct.num_kpts * ct.num_states * P0_BASIS + 1024) * sizeof(OrbitalType));
    RmgCudaError(__FILE__, __LINE__, custat, "cudaMallocHost failed");
    custat = cudaMallocHost((void **)&nv, ct.num_states * P0_BASIS * sizeof(OrbitalType));
    RmgCudaError(__FILE__, __LINE__, custat, "cudaMallocHost failed");
    custat = cudaMallocHost((void **)&ns, ct.num_states * P0_BASIS * sizeof(OrbitalType));
    RmgCudaError(__FILE__, __LINE__, custat, "cudaMallocHost failed");
    custat = cudaMallocHost((void **)&Bns, ct.num_states * P0_BASIS * sizeof(OrbitalType));
    RmgCudaError(__FILE__, __LINE__, custat, "cudaMallocHost failed");
#else
    // Wavefunctions are actually stored here
    rptr = new OrbitalType[ct.num_kpts * ct.num_states * P0_BASIS + 1024];
    nv = new OrbitalType[ct.num_states * P0_BASIS]();
    ns = new OrbitalType[ct.num_states * P0_BASIS]();
    Bns = new OrbitalType[ct.num_states * P0_BASIS]();
#endif
    pct.nv = (double *)nv;
    pct.ns = (double *)ns;
    pct.Bns = (double *)Bns;


    if((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0)) {
        rptr1 = new double[(ct.num_states + 1) * (P0_BASIS + 4) + 1024];
    }

    kpt1 = ct.num_kpts;
    if (verify ("calculation_mode", "Band Structure Only"))
        kpt1 = 1;
    for (kpt = 0; kpt < kpt1; kpt++)
    {

        Kptr[kpt]->set_pool(rptr);
        Kptr[kpt]->nv = nv;
        Kptr[kpt]->ns = ns;
        Kptr[kpt]->Bns = Bns;

        for (st1 = 0; st1 < ct.num_states; st1++)
        {
            Kptr[kpt]->Kstates[st1].kidx = kpt;
            Kptr[kpt]->Kstates[st1].psi = rptr;
            Kptr[kpt]->Kstates[st1].dvhxc = rptr1;
            Kptr[kpt]->Kstates[st1].vxc = vxc;
            Kptr[kpt]->Kstates[st1].vh = vh;
            Kptr[kpt]->Kstates[st1].vnuc = vnuc;
            Kptr[kpt]->Kstates[st1].pbasis =P0_BASIS;
            Kptr[kpt]->Kstates[st1].sbasis = (PX0_GRID + 4) * (PY0_GRID + 4) * (PZ0_GRID + 4);
            Kptr[kpt]->Kstates[st1].istate = st1;
            Kptr[kpt]->Kstates[st1].vel = get_vel();
            rptr +=P0_BASIS;
            rptr1 +=P0_BASIS;
        }
    }


    //Dprintf ("If not an initial run read data from files");
    if (ct.runflag == RESTART)
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

        /* Set initial states to random start */
    
	if (ct.runflag == LCAO_START)
	{
	    lcao_init ();
            LcaoGetPsi(Kptr[0]->Kstates);
	}
	
	else
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


    }                           /* end if */


    //Dprintf ("Initialize the radial potential stuff");
    init_psp ();

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

	printf ("\n\n init: Starting FFTW initialization ...");

    fflush (NULL);

    /*Setup some fftw stuff */
    /*Setup fftw wisdom */
    //init_fftw_wisdom ();

    //Dprintf ("Get memory for fourier transform phase");
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];

        iptr->fftw_phase_sin = new double[sp->nldim * sp->nldim * sp->nldim];
        iptr->fftw_phase_cos = new double[sp->nldim * sp->nldim * sp->nldim];
    }


    /*Do forward transform for each species and store results on the coarse grid */
    init_weight ();
    /*The same for derivative of beta */
    init_derweight ();

    printf ("\n init: FFTW initialization finished, it took %.1f s", my_crtc () - time2);
    fflush (NULL);


    /* Initialize the qfunction stuff */
    init_qfunct ();

    
    /* Update items that change when the ionic coordinates change */
    ReinitIonicPotentials (Kptr, vnuc, rhocore, rhoc);


    // Normalize orbitals if not an initial run
    if (ct.runflag != RESTART) /* Initial run */
    {
        for (kpt = 0; kpt < ct.num_kpts; kpt++)
        {
            for (state = 0; state < ct.num_states; state++)
            {
                Kptr[kpt]->Kstates[state].normalize(Kptr[kpt]->Kstates[state].psi, state);
            }
        }
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
    
    if (ct.runflag == LCAO_START)
	lcao_get_rho(rho);


    ct.rms = 0.0;
    
    //Dprintf ("Generate initial vxc potential and hartree potential");
    pack_vhstod (vh, ct.vh_ext, FPX0_GRID, FPY0_GRID, FPZ0_GRID, ct.boundaryflag);



    //Dprintf ("Condition of run flag is %d", ct.runflag);
    /*If not a restart, get vxc and vh, for restart these quantities should read from the restart file*/
    if (ct.runflag != RESTART)
    {
       	get_vxc (rho, rho_oppo, rhocore, vxc);

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
        Betaxpsi (Kptr[kpt]);
        Kptr[kpt]->mix_betaxpsi(0);
    }

    // If not a restart and diagonalization is requested do a subspace diagonalization otherwise orthogonalize
    if(ct.runflag != RESTART) {

        if (ct.initdiag)
        {
            /*dnmI has to be stup before calling subdiag */
            vtot = new double[FP0_BASIS];

            for (idx = 0; idx < FP0_BASIS; idx++)
                vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];

            /*Generate the Dnm_I */
            get_ddd (vtot);

            /*Release vtot memory */
            delete [] vtot;

            /*Now we cam do subspace diagonalization */
            for(kpt = 0;kpt < ct.num_kpts;kpt++) {
                Subdiag (Kptr[0], vh, vnuc, vxc, ct.subdiag_driver);
            }

        }                           /*end if(ct.initdiag) */
        else
        {

            for (kpt = 0; kpt < ct.num_kpts; kpt++) {
                Kptr[kpt]->orthogonalize(Kptr[kpt]->orbital_storage);
            }
            
        }

    }


#if 0
    /*Setup things for mulliken population analysis */
    if (ct.domilliken)
    {

        int i, it1;
        double t1, t2, scale;
        /* Find milliken radius in number of grid points */
        /* Set the scaling factor for determining the radius of the local grids */
        scale = 1.0;
        if (ct.ibrav == CUBIC_BC)
            scale = 1.1;
        if (ct.ibrav == CUBIC_FC)
            scale = 1.3;

        /*Loop over species */
        for (species = 0; species < ct.num_species; species++)
        {

            /* Get species type */
            sp = &ct.sp[species];

            t1 = 2.0 * scale * sp->mill_radius / ct.hmingrid;
            t1 = modf (t1, &t2);
            it1 = (int) t2;
            if (t1 > 0.5)
                it1++;
            if (!(it1 % 2))
                it1++;
            sp->mill_dim = it1;

            if ((sp->mill_dim >= get_NX_GRID()) || (sp->mill_dim >= get_NY_GRID())
                || (sp->mill_dim >= get_NZ_GRID()))
                error_handler ("Milliken radius exceeds global grid size");

            /*Find total number of states -i.e. we only know that we have states s,p, d etc
             * but we need number of s, px, py, pz, etc*/
            sp->sum_atomic_waves = 0;

            i = 0;
            for (st1 = 0; st1 < sp->num_atomic_waves; st1++)
            {
                sp->sum_atomic_waves += 2 * sp->lstate_atomic_wave[st1] + 1;

                switch (sp->lstate_atomic_wave[st1])
                {
                case 0:
                    strcpy (sp->atomic_wave_symbol[i], "s");
                    i++;
                    break;

                case 1:
                    strcpy (sp->atomic_wave_symbol[i], "px");
                    strcpy (sp->atomic_wave_symbol[i + 1], "pz");
                    strcpy (sp->atomic_wave_symbol[i + 2], "py");
                    i += 3;
                    break;

                case 2:
                    strcpy (sp->atomic_wave_symbol[i], "dxy");
                    strcpy (sp->atomic_wave_symbol[i + 1], "dxz");
                    strcpy (sp->atomic_wave_symbol[i + 2], "dzz");
                    strcpy (sp->atomic_wave_symbol[i + 3], "dyz");
                    strcpy (sp->atomic_wave_symbol[i + 4], "dxx-yy");
                    i += 5;
                    break;

                case 3:
                    strcpy (sp->atomic_wave_symbol[i], "Fxxx");
                    strcpy (sp->atomic_wave_symbol[i + 1], "Fyyy");
                    strcpy (sp->atomic_wave_symbol[i + 2], "Fxyz");
                    strcpy (sp->atomic_wave_symbol[i + 3], "Fzzz");
                    strcpy (sp->atomic_wave_symbol[i + 4], "Fz(xx-yy)");
                    strcpy (sp->atomic_wave_symbol[i + 5], "Fy(zz-xx)");
                    strcpy (sp->atomic_wave_symbol[i + 6], "Fx(yy-zz)");
                    i += 7;
                }               /*end switch */

            }



        }                       /*end loop over species */

    }                           /*end if (ct.domilliken) */
#endif




}                               /* end init */




static void init_alloc_nonloc_mem (void)
{
    int ion;

#if FDIFF_BETA
    pct.weight_derx = new double *[ct.num_ions];
    pct.weight_dery = new double *[ct.num_ions];
    pct.weight_derz = new double *[ct.num_ions];
#endif

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

#if FDIFF_BETA
        pct.weight_derx[ion] = NULL;
        pct.weight_dery[ion] = NULL;
        pct.weight_derz[ion] = NULL;
#endif

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
