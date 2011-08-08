/************************** SVN Revision Information **************************
 **    $Id$    **
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
 *   void init(REAL *vh, REAL *rho, REAL *rhocore, REAL *rhoc, 
 *             STATE *states, REAL *vnuc, REAL *vxc)
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
#include "main.h"


static void init_alloc_nonloc_mem (void);



void init (REAL * vh, REAL * rho, REAL * rho_oppo, REAL * rhocore, REAL * rhoc,
           STATE * states, REAL * vnuc, REAL * vxc)
{

    int kpt, kpt1, kst1, ic, idx, state, ion, st1, it1;
    int species, i;
    SPECIES *sp;
    int flag;
    REAL *rptr = NULL, *vtot, t1, t2, scale, *rho_tot;
    ION *iptr, *iptr0;
    REAL time1, time2, v1, v2, v3, fac;
    STATE *st;

#ifdef SMP
    int thread, offset, stop;
#endif

    time1 = my_crtc ();

    ct.fftw_wisdom_setup = 0;

    ct.total_scf_steps = 0;
    ct.md_steps = 0;

    ct.psi_nbasis = NX_GRID * NY_GRID * NZ_GRID;
    ct.psi_nxgrid = NX_GRID;
    ct.psi_nygrid = NY_GRID;
    ct.psi_nzgrid = NZ_GRID;

    ct.psi_fnbasis = FNX_GRID * FNY_GRID * FNZ_GRID;
    ct.psi_fnxgrid = FNX_GRID;
    ct.psi_fnygrid = FNY_GRID;
    ct.psi_fnzgrid = FNZ_GRID;


    /*Initialize ScaLapack, get context */
    sl_init (&pct.ictxt, ct.num_states);

    /*Set desca, variable used in ScaLapack functions */
    set_desca (pct.desca, &pct.ictxt, ct.num_states);


    /* Allocate storage for non-local projectors */
    pct.newsintR_local = NULL;
    pct.oldsintR_local = NULL;
#if !GAMMA_PT
    pct.newsintI_local = NULL;
    pct.oldsintI_local = NULL;
#endif


    /* Set hartree boundary condition stuff */
    ct.vh_pxgrid = FPX0_GRID;
    ct.vh_pygrid = FPY0_GRID;
    ct.vh_pzgrid = FPZ0_GRID;

    ct.vh_pbasis = ct.vh_pxgrid * ct.vh_pygrid * ct.vh_pzgrid;
    my_malloc (ct.vh_ext, ct.vh_pbasis, REAL);


    /* initialize the lattice basis vectors */
    flag = 0;
    latgen (&ct.ibrav, ct.celldm, ct.a0, ct.a1, ct.a2, &ct.omega, &flag);


    /* initialize the reciprocal lattice vectors */
    recips ();

    /* Initialize some k-point stuff */
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        v1 = twoPI * ct.kp[kpt].kpt[0] / ct.xside;
        v2 = twoPI * ct.kp[kpt].kpt[1] / ct.yside;
        v3 = twoPI * ct.kp[kpt].kpt[2] / ct.zside;

        ct.kp[kpt].kvec[0] = v1;
        ct.kp[kpt].kvec[1] = v2;
        ct.kp[kpt].kvec[2] = v3;

        ct.kp[kpt].kmag = v1 * v1 + v2 * v2 + v3 * v3;

        ct.kp[kpt].kstate = &states[kpt * ct.num_states];
        ct.kp[kpt].kidx = kpt;
    }


    /* Get the crystal or cartesian coordinates of the ions */
    init_pos ();


    ct.hmaxgrid = ct.xside * ct.hxgrid;
    if (ct.yside * ct.hygrid > ct.hmaxgrid)
        ct.hmaxgrid = ct.yside * ct.hygrid;
    if (ct.zside * ct.hzgrid > ct.hmaxgrid)
        ct.hmaxgrid = ct.zside * ct.hzgrid;

    ct.hmingrid = ct.xside * ct.hxgrid;
    if (ct.yside * ct.hygrid < ct.hmingrid)
        ct.hmingrid = ct.yside * ct.hygrid;
    if (ct.zside * ct.hzgrid < ct.hmingrid)
        ct.hmingrid = ct.zside * ct.hzgrid;

    ct.anisotropy = ct.hmaxgrid / ct.hmingrid;

    if (ct.anisotropy > 1.1)
    {
        if (pct.imgpe == 0)
        {
            printf ("ct.hxgrid = %7.5f\n", ct.hxgrid);
            printf ("ct.hygrid = %7.5f\n", ct.hygrid);
            printf ("ct.hzgrid = %7.5f\n", ct.hzgrid);
        }
        error_handler ("Anisotropy too large");
    }


    /* Set discretization array */
    ct.xcstart = 0.0;
    ct.ycstart = 0.0;
    ct.zcstart = 0.0;


    /* Some multigrid parameters */
    ct.poi_parm.sb_step = 1.0;
    ct.eig_parm.sb_step = 1.0;

#if MPI
    /* Set state pointers and initialize state data */
#if GAMMA_PT
    /* Wavefunctions are actually stored here */
    my_malloc (rptr, (ct.num_states + 1) * (P0_BASIS + 4) + 1024, REAL);
#else
    /* Wavefunctions are actually stored here */
    if (verify ("calculation_mode", "Band Structure Only"))
    {
        my_malloc (rptr, 2 * (ct.num_states + 1) * (P0_BASIS + 4) + 1024, REAL);
    }
    else
    {
        my_malloc (rptr, ct.num_kpts * 2 * (ct.num_states + 1) * (P0_BASIS + 4) + 1024, REAL);
    }
#endif

#endif
    kpt1 = ct.num_kpts;
    if (verify ("calculation_mode", "Band Structure Only"))
        kpt1 = 1;
    kst1 = 0;
    for (kpt = 0; kpt < kpt1; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {
            states[kst1].kidx = kpt;
            states[kst1].psiR = rptr;
            states[kst1].psiI = rptr + P0_BASIS;
            states[kst1].hxgrid = ct.hxgrid;
            states[kst1].hygrid = ct.hygrid;
            states[kst1].hzgrid = ct.hzgrid;
            states[kst1].dimx = PX0_GRID;
            states[kst1].dimy = PY0_GRID;
            states[kst1].dimz = PZ0_GRID;
            states[kst1].vxc = vxc;
            states[kst1].vh = vh;
            states[kst1].vnuc = vnuc;
            states[kst1].pbasis = P0_BASIS;
            states[kst1].sbasis = (PX0_GRID + 2) * (PY0_GRID + 2) * (PZ0_GRID + 2);
            states[kst1].istate = st1;
            states[kst1].vel = ct.vel;
#if MPI
#if GAMMA_PT
            rptr += P0_BASIS;
#else
            rptr += 2 * P0_BASIS;
#endif
#endif
            kst1++;
        }
    }


    Dprintf ("If not an initial run read data from files");
    if (ct.runflag == 1)
    {
        read_data (ct.infile, vh, rho, rho_oppo, vxc, states);
    }
    else if ((ct.runflag == 0) || (ct.runflag == 2))    /* Initial run */
    {
        /* Set the initial hartree potential to a constant */
        for (idx = 0; idx < FP0_BASIS; idx++)
            vh[idx] = 0.0;

        /* Set initial states to random start */
        for (kpt = 0; kpt < ct.num_kpts; kpt++)
            init_wf (&states[kpt * ct.num_states]);


        /* Set initial ionic coordinates to the current ones. */
        for (ion = 0; ion < ct.num_ions; ion++)
        {
            for (ic = 0; ic < 3; ic++)
            {
                ct.ions[ion].icrds[ic] = ct.ions[ion].crds[ic];
                ct.ions[ion].ixtal[ic] = ct.ions[ion].xtal[ic];
            }
        }


        /* set the initial force pointer array up */
        ct.fpt[0] = 0;
        ct.fpt[1] = 1;
        ct.fpt[2] = 2;
        ct.fpt[3] = 3;


        /*Initialize velocities */
        for (ion = 0; ion < ct.num_ions; ion++)
            for (ic = 0; ic < 3; ic++)
                ct.ions[ion].velocity[ic] = 0.0;

    }                           /* end if */


    Dprintf ("Initialize the radial potential stuff");
    init_psp ();

    /* Initialize symmetry stuff */
#if !GAMMA_PT
    init_sym ();
#endif


    Dprintf ("Allocate memory for arrays related to nonlocal PP");
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
    init_fftw_wisdom ();

    Dprintf ("Get memory for fourier transform phase");
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];

        my_malloc (iptr->fftw_phase_sin, sp->nldim * sp->nldim * sp->nldim, REAL);
        my_malloc (iptr->fftw_phase_cos, sp->nldim * sp->nldim * sp->nldim, REAL);
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
    reinit_ionic_pp (states, vnuc, rhocore, rhoc);


    if ((ct.runflag == 0) || (ct.runflag == 2)) /* Initial run */
    {
        for (kpt = 0; kpt < ct.num_kpts; kpt++)
        {
            for (state = 0; state < ct.num_states; state++)
            {
                st = &states[kpt * ct.num_states + state];
                norm_psi1 (st, state, kpt);
            }
        }
    }
        
    mix_betaxpsi(0);

    Dprintf ("Set the initial density to be equal to the compensating charges");
    if (!ct.runflag)
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


#ifdef SMP
    offset = 0;
    for (thread = 0; thread < ct.num_threads; thread++)
    {
        stop = thread_control[thread].numpt;
        rptr = thread_control[thread].rho;
        QMD_scopy (stop, &rho[offset], 1, rptr, 1);
        offset += stop;
    }
#endif


    Dprintf ("Generate initial vxc potential and hartree potential");
    pack_vhstod (vh, ct.vh_ext, FPX0_GRID, FPY0_GRID, FPZ0_GRID);



    Dprintf ("Condition of run flag is %d", ct.runflag);
    if (ct.runflag != 1)
    {
       	get_vxc (rho, rho_oppo, rhocore, vxc);
       	Dprintf ("get vxc completed");

	if (ct.spin_flag)
	{
    		my_malloc (rho_tot,  FP0_BASIS, REAL);
  	    	for (idx = 0; idx < FP0_BASIS; idx++)
            		rho_tot[idx] = rho[idx] + rho_oppo[idx];

        	get_vh (rho_tot, rhoc, vh, 30, ct.poi_parm.levels);
    		my_free (rho_tot);
	}
	else
        	get_vh (rho, rhoc, vh, 30, ct.poi_parm.levels);

    }

    Dprintf ("If diagonalization is requested do a subspace diagonalization");
    if (ct.initdiag)
    {
        /*dnmI has to be stup before calling subdiag */
        my_malloc (vtot, FP0_BASIS, REAL);


        for (idx = 0; idx < FP0_BASIS; idx++)
            vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];

        /*Generate the Dnm_I */
        get_ddd (vtot);

        /*Release vtot memory */
        my_free (vtot);

        /*Now we cam do subspace diagonalization */
#if GAMMA_PT
        subdiag_gamma (states, vh, vnuc, vxc);
#else
        subdiag_nongamma (states, vh, vnuc, vxc);
#endif

    }                           /*end if(ct.initdiag) */



    /*Setup things for mulliken population analysis */
    if (ct.domilliken)
    {

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

            if ((sp->mill_dim >= ct.psi_nxgrid) || (sp->mill_dim >= ct.psi_nygrid)
                || (sp->mill_dim >= ct.psi_nzgrid))
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


    rmg_timings (INIT_TIME, (my_crtc () - time1));


}                               /* end init */




static void init_alloc_nonloc_mem (void)
{
    int ion;

    my_malloc (pct.weight, ct.num_ions, REAL *);

#if FDIFF_BETA
    my_malloc (pct.weight_derx, ct.num_ions, REAL *);
    my_malloc (pct.weight_dery, ct.num_ions, REAL *);
    my_malloc (pct.weight_derz, ct.num_ions, REAL *);
#endif

    my_malloc (pct.nlindex, ct.num_ions, int *);
    my_malloc (pct.Qindex, ct.num_ions, int *);
    my_malloc (pct.idxflag, ct.num_ions, int *);
    my_malloc (pct.Qdvec, ct.num_ions, int *);


    my_malloc (pct.idxptrlen, ct.num_ions, int);
    my_malloc (pct.Qidxptrlen, ct.num_ions, int);
    my_malloc (pct.lptrlen, ct.num_ions, int);


    my_malloc (pct.phaseptr, ct.num_ions, REAL *);
    my_malloc (pct.augfunc, ct.num_ions, REAL *);
    my_malloc (pct.dnmI, ct.num_ions, REAL *);
    my_malloc (pct.qqq, ct.num_ions, REAL *);

    //my_malloc(pct.nonloc_ions_list, ct.num_ions, int);
    //my_malloc(pct.q_ions_list, ct.num_ions, int);


    /*Initialize pointer arrays to NULL */
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        pct.weight[ion] = NULL;

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
