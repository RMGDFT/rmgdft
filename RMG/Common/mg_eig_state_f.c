/************************** SVN Revision Information **************************
 **    $Id: mg_eig_state.c 1908 2013-02-25 12:51:00Z ebriggs $    **
******************************************************************************/

/****f* QMD-MGDFT/mg_eig_state.c *****
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
 *   void mg_eig_state_f(STATE *sp, int tid, REAL *vtot)
 *   mixed precision gamma point only  Multigrid solver for a single electronic orbital.
 *
 *   Algorithm:
 *
 *     1.  Eigenvalues are computed once per cycle using the Mehrstallen
 *         operator.
 *
 *     2.  Presmoothing of wavefunctions on the finest grid level using
 *         the Mehrstallen operator.
 *
 *     3.  The defect from the last Mehrstallen smoothing is then used
 *         as the right hand side for a richardson defect correction
 *         scheme using a 7-point discretization of the del-squared operator.
 *         This is multigridded.
 *
 *     4.  Mehrstallen post smoothing.
 * INPUTS
 *   sp: point to orbital structure STATE (see main.h)
 *   tid: thread ID
 *   vtot: total potential = vxc +vnuc +vh
 * OUTPUT
 *   wave functions and eigenvalues are updated and stored in STATE
 * PARENTS
 *   scf.c threads.c
 * CHILDREN
 *   gather_psi.c app_nl.c pack_ptos.c app_cilr.c genvpsi.c app_cir.c global_sums.c
 *   trade_images.c app_smooth.c mgrid_solv.c pack_stop_axpy.c scatter_psi.c
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


#if HYBRID_MODEL
    #include "hybrid.h"
    #include <pthread.h>
static pthread_mutex_t vtot_sync_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif


extern STATE *states;

void mg_eig_state_f (STATE * sp, int tid, REAL * vtot_psi)
{

    int idx, cycles, ntid;
    int nits, pbasis, sbasis;
    REAL eig, diag, t1, t2, t3, t4;
    REAL *work1, *nv, *ns, *res2;
    REAL *tmp_psi, *res, *saved_psi;
    REAL *nvtot_psi;
    rmg_float_t *tmp_psi_f, *work2_f, *res_f, *res2_f, *sg_psi_f;
    int eig_pre[6] = { 0, 3, 6, 2, 2, 2 };
    int eig_post[6] = { 0, 3, 6, 2, 2, 2 };
    int ione = 1;
    int dimx, dimy, dimz, levels, potential_acceleration;
    REAL hxgrid, hygrid, hzgrid, sb_step;
    REAL tarr[8];
    REAL time1;
    rmg_float_t *sg_twovpsi_f, *work1_f;

    nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;
    dimx = sp->dimx;
    dimy = sp->dimy;
    dimz = sp->dimz;
    hxgrid = sp->hxgrid;
    hygrid = sp->hygrid;
    hzgrid = sp->hzgrid;
    levels = ct.eig_parm.levels;
    sb_step = 1.0;
    pbasis = sp->pbasis;
    sbasis = sp->sbasis;

    /* Grab some memory */
    my_malloc (sg_psi_f, sbasis, rmg_float_t);
    my_malloc (res, sbasis, REAL);
#if GPU_ENABLED
#if HYBRID_MODEL
    ntid = get_thread_tid();
#else
    ntid = 0;
#endif
    if(ntid == -1) {         // Normal codepath with no threads
        work2_f = (rmg_float_t *)&ct.gpu_host_temp2[0];
        sg_twovpsi_f = (rmg_float_t *)&ct.gpu_host_temp1[0];
        work1 = &ct.gpu_host_temp4[0];
    }
    else {                  // Threaded codepath for hybrid mode since each thread needs it's own copy
        work2_f = (rmg_float_t *)&ct.gpu_host_temp2[ntid * 4 * sbasis];
        sg_twovpsi_f = (rmg_float_t *)&ct.gpu_host_temp1[ntid * sbasis];
        work1 = &ct.gpu_host_temp4[ntid * sbasis];
    }

#else
    my_malloc (work2_f, 4 * sbasis, rmg_float_t);
    my_malloc (sg_twovpsi_f, sbasis, rmg_float_t);
    my_malloc (work1, sbasis, REAL);
#endif
    my_malloc (ns, sbasis, REAL);
    my_malloc (nv, sbasis, REAL);
    my_malloc (res2, sbasis, REAL);
    my_malloc (saved_psi, sbasis, REAL);
    my_malloc (nvtot_psi, sbasis, REAL);
    my_malloc (work1_f, 4*sbasis, rmg_float_t);
    my_malloc (tmp_psi_f, sbasis, rmg_float_t);
    my_malloc (res_f, sbasis, rmg_float_t);
    my_malloc (res2_f, sbasis, rmg_float_t);

    tmp_psi = sp->psiR;

    potential_acceleration = ((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0));
    if(potential_acceleration) {
        for(idx = 0;idx <pct.P0_BASIS;idx++) {
            nvtot_psi[idx] = vtot_psi[idx];
            saved_psi[idx] = tmp_psi[idx];
        }

    }

#if MD_TIMERS
    time1 = my_crtc ();
#endif

    if(ct.eig_parm.mucycles > 1)
        mix_betaxpsi1(sp);

    /* Get the non-local operator and S acting on psi (nv and ns, respectively) */
    app_nls (tmp_psi, NULL, nv, NULL, ns, NULL, pct.oldsintR_local, NULL, sp->istate, sp->kidx);

    // Copy double precision psi into single precison array
    for(idx = 0;idx < pbasis;idx++)
        tmp_psi_f[idx] = (rmg_float_t)tmp_psi[idx];

    // Copy double precision ns into temp single precision array */
    for(idx = 0;idx < pbasis;idx++)
        work1_f[idx] = (rmg_float_t)ns[idx];
  
#if MD_TIMERS
    rmg_timings (MG_EIG_NLS_TIME, (my_crtc () - time1));
#endif

#if MD_TIMERS
    time1 = my_crtc ();
#endif
    /*Apply double precision Mehrstellen right hand operator to ns (stored temporarily in work1_f) and save in res2 */
    app_cir_driver_f (work1_f, res2_f, dimx, dimy, dimz, ct.kohn_sham_fd_order);

#if GPU_ENABLED
    //cudaDeviceSynchronize();
#endif

#if MD_TIMERS
    rmg_timings (MG_EIG_APPCIR_TIME, (my_crtc () - time1));
#endif


    /* Smoothing cycles */
    for (cycles = 0; cycles <= nits; cycles++)
    {

#if MD_TIMERS
        time1 = my_crtc ();
#endif

        //scf_barrier_wait();
        /* Apply Mehrstellen left hand operators */
        diag = app_cil_driver_f (tmp_psi_f, work2_f, dimx, dimy, dimz, hxgrid, hygrid, hzgrid, ct.kohn_sham_fd_order);

#if MD_TIMERS
        rmg_timings (MG_EIG_APPCIL_TIME, (my_crtc () - time1));
#endif

        diag = -1.0 / diag;

#if MD_TIMERS
        time1 = my_crtc ();
#endif

        // Copy saved application to ns to res
        QMD_scopy(pbasis, res2_f, 1, res_f, 1);

#if MD_TIMERS
        time1 = my_crtc ();
#endif

        if(potential_acceleration) {
            /* Generate 2 * V * psi */
            genvpsi_f (tmp_psi_f, sg_twovpsi_f, nvtot_psi, nv, NULL, 0.0, dimx, dimy, dimz);
        }
        else { 
            genvpsi_f (tmp_psi_f, sg_twovpsi_f, vtot_psi, nv, NULL, 0.0, dimx, dimy, dimz);
        }

#if MD_TIMERS
        rmg_timings (MG_EIG_GENVPSI_TIME, (my_crtc () - time1));
#endif

#if MD_TIMERS
        time1 = my_crtc ();
#endif
#if GPU_ENABLED
        cudaDeviceSynchronize();
#endif

        //scf_barrier_wait();
        /* B operating on 2*V*psi stored in work1 */
        app_cir_driver_f (sg_twovpsi_f, work1_f, dimx, dimy, dimz, ct.kohn_sham_fd_order);

#if GPU_ENABLED
        //cudaDeviceSynchronize();
#endif

#if MD_TIMERS
        rmg_timings (MG_EIG_APPCIR_TIME, (my_crtc () - time1));
#endif

        t1 = -ONE;
        QMD_saxpy (pbasis, t1, work2_f, ione, work1_f, ione);


#if MD_TIMERS
        time1 = my_crtc ();
#endif
        /* If this is the first time through compute the eigenvalue */
        if ((cycles == 0) || (potential_acceleration != 0)) 
        {

            eig = 0.0;
            t2 = 0.0;

            for (idx = 0; idx < pbasis; idx++)
            {

                t2 += (rmg_double_t)tmp_psi_f[idx] * (rmg_double_t)res_f[idx];
                eig += (rmg_double_t)tmp_psi_f[idx] * (rmg_double_t)work1_f[idx];

            }

            tarr[0] = t2;
            tarr[1] = eig;
            idx = 2;
            global_sums (tarr, &idx, pct.grid_comm);


            /*If diagonalization is done every step, do not calculate eigenvalues, use those
             * from diagonalization, except for the first step, since at that time eigenvalues 
	     * are not defined yet*/
            if ((ct.diag == 1) && (potential_acceleration == 0) && (ct.scf_steps < ct.end_diag))
            {
                if (ct.scf_steps == 0)
                {
                    eig = tarr[1] / (TWO * tarr[0]);
                    sp->eig[0] = eig;
                    sp->oldeig[0] = eig;
                }
		else
                    eig = sp->eig[0];
	    }
            else
            {
                eig = tarr[1] / (TWO * tarr[0]);
                sp->eig[0] = eig;
                if(ct.scf_steps == 0) {
                    sp->oldeig[0] = eig;
                }
            }
            
            if(potential_acceleration) {
                t1 = eig;
                eig = 0.3 * eig + 0.7 * sp->oldeig[0];
                sp->oldeig[0] = t1;
            }

        }

#if MD_TIMERS
        rmg_timings (MG_EIG_EIGVALUE_TIME, (my_crtc () - time1));
#endif
//dprintf("EIG0 = %d  %18.12f", sp->istate,27.2*sp->eig[0]);

        /* Now either smooth the wavefunction or do a multigrid cycle */
        if (cycles == ct.eig_parm.gl_pre)
        {


            t1 = TWO * eig;
            for (idx = 0; idx <pct.P0_BASIS; idx++)
            {

                res_f[idx] = t1 * res_f[idx] - work1_f[idx];

            }

            /* Pack the residual data into multigrid array */
            pack_ptos_f (sg_psi_f, res_f, dimx, dimy, dimz);

#if MD_TIMERS
            time1 = my_crtc ();
#endif
            trade_images_f (sg_psi_f, dimx, dimy, dimz, pct.neighbors);

#if MD_TIMERS
            rmg_timings (MG_EIG_TRADE_TIME, (my_crtc () - time1));
#endif

            /* Smooth it once and store the smoothed residual in work1 */
            t1 = 145.0;

#if MD_TIMERS
            time1 = my_crtc ();
#endif
            app_smooth_f (sg_psi_f, work1_f, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);
#if MD_TIMERS
            rmg_timings (MG_EIG_APPSMOOTH_TIME, (my_crtc () - time1));
#endif

#if MD_TIMERS
            time1 = my_crtc ();
#endif

            if(potential_acceleration) {
                t1 = eig - states[0].eig[0];
                t1 = -t1*t1 / 10.0;
            }
            else {
                t1 = 0.0;
            }

            /* Do multigrid step with solution returned in sg_twovpsi */
            mgrid_solv_f (sg_twovpsi_f, work1_f, work2_f,
                        dimx, dimy, dimz, hxgrid,
                        hygrid, hzgrid, 0, pct.neighbors, levels, eig_pre, eig_post, 1, sb_step, t1,
                        NX_GRID, NY_GRID, NZ_GRID,
                        pct.PX_OFFSET, pct.PY_OFFSET, pct.PZ_OFFSET,
                        pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);

#if MD_TIMERS
            rmg_timings (MG_EIG_MGRIDSOLV_TIME, (my_crtc () - time1));
#endif


            /* The correction is in a smoothing grid so we use this
             * routine to update the orbital which is stored in a physical grid.
             */
#if MD_TIMERS
            time1 = my_crtc ();
#endif

            t1 = -1.0;
            pack_stop_axpy_f (sg_twovpsi_f, tmp_psi_f, t1, dimx, dimy, dimz);

#if MD_TIMERS
            rmg_timings (MG_EIG_PACK_TIME, (my_crtc () - time1));
#endif

        }
        else
        {


            t1 = TWO * eig;
            t2 = ZERO;
            t4 = ct.eig_parm.gl_step * diag;

            for (idx = 0; idx <pct.P0_BASIS; idx++)
            {

                t3 = t1 * (rmg_double_t)res_f[idx] - (rmg_double_t)work1_f[idx];
                t2 += t3 * t3;
                tmp_psi_f[idx] += t4 * t3;

            }

            if (cycles == 0)
            {

                t2 = real_sum_all (t2, pct.grid_comm);
                t1 = (REAL) (ct.psi_nbasis);
                sp->res = ct.hmaxgrid * ct.hmaxgrid * sqrt (t2 / t1) * 0.25;

            }

        }

    }                           /* end for */


    if(potential_acceleration) {

        // Save potential used for this orbital and update potential for future orbitals
        for(idx = 0;idx <pct.P0_BASIS;idx++) {
            sp->dvhxc[idx] = nvtot_psi[idx];
        }


        if(ct.potential_acceleration_constant_step > 0.0) {

            t1 = 1.8 * ct.potential_acceleration_constant_step;
            if(sp->occupation[0] < 0.5) t1 = 0.0;

#if HYBRID_MODEL
            pthread_mutex_lock(&vtot_sync_mutex);
#endif
            for(idx = 0;idx <pct.P0_BASIS;idx++) {
               vtot_psi[idx] = vtot_psi[idx] + t1 * PI * sp->occupation[0] * tmp_psi_f[idx] * (tmp_psi_f[idx] - saved_psi[idx]);
            }
#if HYBRID_MODEL
            pthread_mutex_unlock(&vtot_sync_mutex);
#endif

        }

        if(ct.potential_acceleration_poisson_step > 0.0) {

            // construct delta_rho
            for(idx = 0;idx <pct.P0_BASIS;idx++) {
                res_f[idx] = -4.0 * PI * sp->occupation[0] *
                           (tmp_psi_f[idx] - saved_psi[idx]) * (2.0*saved_psi[idx] + (tmp_psi_f[idx] - saved_psi[idx]));
            }

            // zero out solution vector
            for(idx = 0;idx <pct.P0_BASIS;idx++) {
                sg_twovpsi_f[idx] = 0.0;
            }

            /* Pack delta_rho into multigrid array */
            pack_ptos_f (sg_psi_f, res_f, dimx, dimy, dimz);
            trade_images_f (sg_psi_f, dimx, dimy, dimz, pct.neighbors);
            /* Smooth it once and store the smoothed charge in res */
            app_smooth1_f (sg_psi_f, res_f, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);

            // neutralize cell with a constant background charge
            t2 = 0.0;
            for(idx = 0;idx <pct.P0_BASIS;idx++) {
                t2 += res_f[idx];
            }
            t2 = real_sum_all(t2, pct.grid_comm) / (NX_GRID * NY_GRID * NZ_GRID);
            for(idx = 0;idx <pct.P0_BASIS;idx++) {
                res_f[idx] -= t2;
            }

            /* Do multigrid step with solution returned in sg_twovpsi */
            eig_pre[0] = 2;
            eig_post[0] = 2;
            levels=1;
            mgrid_solv_f (sg_twovpsi_f, res_f, work2_f,
                        dimx, dimy, dimz, hxgrid,
                        hygrid, hzgrid, 0, pct.neighbors, levels, eig_pre, eig_post, 1, 1.0, 0.0,
                        NX_GRID, NY_GRID, NZ_GRID,
                        pct.PX_OFFSET, pct.PY_OFFSET, pct.PZ_OFFSET,
                        pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);

            for(idx = 0;idx <pct.P0_BASIS;idx++) {
                res_f[idx] = 0.0;
            }
            pack_stop_axpy_f (sg_twovpsi_f, res_f, 1.0, dimx, dimy, dimz);
            t1 = ct.potential_acceleration_poisson_step;
            if(sp->occupation[0] < 0.5) t1 = 0.0;
#if HYBRID_MODEL
            pthread_mutex_lock(&vtot_sync_mutex);
#endif
            for(idx = 0;idx <pct.P0_BASIS;idx++) {
               vtot_psi[idx] = vtot_psi[idx] + t1 * res_f[idx];
            }
#if HYBRID_MODEL
            pthread_mutex_unlock(&vtot_sync_mutex);
#endif
        }


    } // end if

    // Copy single precision orbital back to double precision
    for(idx = 0;idx < pbasis;idx++) {

        tmp_psi[idx] = (rmg_double_t)tmp_psi_f[idx];

    }

//dprintf("EIG1 = %d  %18.12f", sp->istate,27.2*sp->eig[0]);
    /* Release our memory */
    my_free (res_f);
    my_free (res2_f);
    my_free (tmp_psi_f);
    my_free (work1_f);
    my_free (nvtot_psi);
    my_free (saved_psi);
    my_free (res2);
    my_free (nv);
    my_free (ns);
#if !GPU_ENABLED
    my_free (work2_f);
    my_free (sg_twovpsi_f);
    my_free (work1);
#endif
    my_free (res);
    my_free (sg_psi_f);


}                               /* end mg_eig_state_f */

