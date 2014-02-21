/************************** SVN Revision Information **************************
 **    $Id$    **
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
 *   void mg_eig_state(STATE *sp, int tid, rmg_double_t *vtot)
 *   there are two routines, one for real (Gamma point only) and the other for complex
 *   Multigrid solver for a single electronic orbital.
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
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"


#if HYBRID_MODEL
    #include "hybrid.h"
    #include <pthread.h>
static pthread_mutex_t vtot_sync_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif


extern STATE *states;

#if GAMMA_PT
void mg_eig_state (STATE * sp, int tid, rmg_double_t * vtot_psi)
{

    int idx, cycles, P0_BASIS;
    int nits, pbasis, sbasis;
    rmg_double_t eig, diag, t1, t2, t3, t4;
    rmg_double_t *work1, *work2, *nv, *ns, *res2;
    rmg_double_t *tmp_psi, *res, *sg_psi, *sg_twovpsi, *saved_psi;
    rmg_double_t *nvtot_psi;
    int eig_pre[6] = { 0, 3, 6, 2, 2, 2 };
    int eig_post[6] = { 0, 3, 6, 2, 2, 2 };
    int ione = 1;
    int dimx, dimy, dimz, levels, potential_acceleration;

    rmg_double_t hxgrid, hygrid, hzgrid, sb_step;
    rmg_double_t tarr[8];
    rmg_double_t time1;

    P0_BASIS = get_P0_BASIS();

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
    my_malloc (sg_psi, sbasis, rmg_double_t);
    my_malloc (res, sbasis, rmg_double_t);

    my_malloc (work2, 4 * sbasis, rmg_double_t);
    my_malloc (sg_twovpsi, 4 * sbasis, rmg_double_t);
    my_malloc (work1, 4 * sbasis, rmg_double_t);

#if !BATCH_NLS
    my_malloc (ns, sbasis, rmg_double_t);
    my_malloc (nv, sbasis, rmg_double_t);
#endif
    my_malloc (res2, sbasis, rmg_double_t);
    my_malloc (saved_psi, sbasis, rmg_double_t);
    my_malloc (nvtot_psi, sbasis, rmg_double_t);

    tmp_psi = sp->psiR;

    potential_acceleration = ((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0));
    if(potential_acceleration) {
        for(idx = 0;idx <P0_BASIS;idx++) {
            nvtot_psi[idx] = vtot_psi[idx];
            saved_psi[idx] = tmp_psi[idx];
        }

    }

#if MD_TIMERS
    time1 = my_crtc ();
#endif

    if(ct.eig_parm.mucycles > 1)
        mix_betaxpsi1(sp);

#if !BATCH_NLS
    app_nls (tmp_psi, NULL, nv, NULL, ns, NULL, pct.oldsintR_local, NULL, sp->istate, sp->kidx);

#else
    nv = &pct.nv[sp->istate * P0_BASIS];
    ns = &pct.ns[sp->istate * P0_BASIS];
#endif


#if MD_TIMERS
    rmg_timings (MG_EIG_NLS_TIME, (my_crtc () - time1));
#endif

#if MD_TIMERS
    time1 = my_crtc ();
#endif
    /*Apply Mehrstellen right hand operator to ns and save in res */
    app_cir_driver (ns, res2, dimx, dimy, dimz, ct.kohn_sham_fd_order);


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
        diag = app_cil_driver (tmp_psi, work2, dimx, dimy, dimz, hxgrid, hygrid, hzgrid, ct.kohn_sham_fd_order);

#if MD_TIMERS
        rmg_timings (MG_EIG_APPCIL_TIME, (my_crtc () - time1));
#endif

        diag = -1.0 / diag;

#if MD_TIMERS
        time1 = my_crtc ();
#endif

        // Copy saved application to ns to res
        QMD_dcopy(pbasis, res2, 1, res, 1);

#if MD_TIMERS
        time1 = my_crtc ();
#endif

        if(potential_acceleration) {
            /* Generate 2 * V * psi */
            genvpsi (tmp_psi, sg_twovpsi, nvtot_psi, nv, NULL, 0.0, dimx, dimy, dimz);
        }
        else { 
            genvpsi (tmp_psi, sg_twovpsi, vtot_psi, nv, NULL, 0.0, dimx, dimy, dimz);
        }

#if MD_TIMERS
        rmg_timings (MG_EIG_GENVPSI_TIME, (my_crtc () - time1));
#endif

#if MD_TIMERS
        time1 = my_crtc ();
#endif

        /* B operating on 2*V*psi stored in work1 */
        app_cir_driver (sg_twovpsi, work1, dimx, dimy, dimz, ct.kohn_sham_fd_order);
        for(idx = 0; idx < dimx * dimy * dimz; idx++) work1[idx] += TWO * nv[idx];

#if MD_TIMERS
        rmg_timings (MG_EIG_APPCIR_TIME, (my_crtc () - time1));
#endif

        t1 = -ONE;
        QMD_daxpy (pbasis, t1, work2, ione, work1, ione);


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

                t2 += tmp_psi[idx] * res[idx];
                eig += tmp_psi[idx] * work1[idx];

            }

            tarr[0] = t2;
            tarr[1] = eig;
            idx = 2;


            /*If diagonalization is done every step, do not calculate eigenvalues, use those
             * from diagonalization, except for the first step, since at that time eigenvalues 
             * are not defined yet*/
            if ((ct.diag == 1) && (potential_acceleration == 0) && (ct.scf_steps < ct.end_diag))
            {
                if (ct.scf_steps == 0)
                {
                    global_sums (tarr, &idx, pct.grid_comm);
                    eig = tarr[1] / (TWO * tarr[0]);
                    sp->eig[0] = eig;
                    sp->oldeig[0] = eig;
                }
                else
                    eig = sp->eig[0];
            }
            else
            {
                global_sums (tarr, &idx, pct.grid_comm);
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

        /* Now either smooth the wavefunction or do a multigrid cycle */
        if (cycles == ct.eig_parm.gl_pre)
        {


            t1 = TWO * eig;
            for (idx = 0; idx <P0_BASIS; idx++)
            {

                res[idx] = t1 * res[idx] - work1[idx];

            }

            /* Pack the residual data into multigrid array */
            pack_ptos (sg_psi, res, dimx, dimy, dimz);

#if MD_TIMERS
            time1 = my_crtc ();
#endif
            trade_images (sg_psi, dimx, dimy, dimz, get_neighbors(), FULL_FD);

#if MD_TIMERS
            rmg_timings (MG_EIG_TRADE_TIME, (my_crtc () - time1));
#endif

            /* Smooth it once and store the smoothed residual in work1 */
            t1 = 145.0;

#if MD_TIMERS
            time1 = my_crtc ();
#endif
            app_smooth (sg_psi, work1, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID());
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
            mgrid_solv (sg_twovpsi, work1, work2,
                        dimx, dimy, dimz, hxgrid,
                        hygrid, hzgrid, 0, get_neighbors(), levels, eig_pre, eig_post, 1, sb_step, t1,
                        get_NX_GRID(), get_NY_GRID(), get_NZ_GRID(),
                        get_PX_OFFSET(), get_PY_OFFSET(), get_PZ_OFFSET(),
                        get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), ct.boundaryflag);

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
            pack_stop_axpy (sg_twovpsi, tmp_psi, t1, dimx, dimy, dimz);

#if MD_TIMERS
            rmg_timings (MG_EIG_PACK_TIME, (my_crtc () - time1));
#endif

        }
        else
        {


            t1 = TWO * eig;
            t2 = ZERO;
            t4 = ct.eig_parm.gl_step * diag;

            for (idx = 0; idx <P0_BASIS; idx++)
            {

                t3 = t1 * res[idx] - work1[idx];
                t2 += t3 * t3;
                tmp_psi[idx] += t4 * t3;

            }

            if (cycles == 0)
            {

                t2 = real_sum_all (t2, pct.grid_comm);
                t1 = (rmg_double_t) (ct.psi_nbasis);
                sp->res = ct.hmaxgrid * ct.hmaxgrid * sqrt (t2 / t1) * 0.25;

            }

        }

    }                           /* end for */


    if(potential_acceleration) {

        // Save potential used for this orbital and update potential for future orbitals
        for(idx = 0;idx <P0_BASIS;idx++) {
            sp->dvhxc[idx] = nvtot_psi[idx];
        }


        if(ct.potential_acceleration_constant_step > 0.0) {

            t1 = 1.8 * ct.potential_acceleration_constant_step;
            if(sp->occupation[0] < 0.5) t1 = 0.0;

#if HYBRID_MODEL
            pthread_mutex_lock(&vtot_sync_mutex);
#endif
            for(idx = 0;idx <P0_BASIS;idx++) {
               vtot_psi[idx] = vtot_psi[idx] + t1 * PI * sp->occupation[0] * tmp_psi[idx] * (tmp_psi[idx] - saved_psi[idx]);
            }
#if HYBRID_MODEL
            pthread_mutex_unlock(&vtot_sync_mutex);
#endif

        }

        if(ct.potential_acceleration_poisson_step > 0.0) {

            // construct delta_rho
            for(idx = 0;idx <P0_BASIS;idx++) {
                res[idx] = -4.0 * PI * sp->occupation[0] *
                           (tmp_psi[idx] - saved_psi[idx]) * (2.0*saved_psi[idx] + (tmp_psi[idx] - saved_psi[idx]));
            }

            // zero out solution vector
            for(idx = 0;idx <P0_BASIS;idx++) {
                sg_twovpsi[idx] = 0.0;
            }

            /* Pack delta_rho into multigrid array */
            pack_ptos (sg_psi, res, dimx, dimy, dimz);
            trade_images (sg_psi, dimx, dimy, dimz, get_neighbors(), FULL_FD);
            /* Smooth it once and store the smoothed charge in res */
            app_smooth1 (sg_psi, res, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID());

            // neutralize cell with a constant background charge
            t2 = 0.0;
            for(idx = 0;idx <P0_BASIS;idx++) {
                t2 += res[idx];
            }
            t2 = real_sum_all(t2, pct.grid_comm) / (get_NX_GRID() * get_NY_GRID() * get_NZ_GRID());
            for(idx = 0;idx <P0_BASIS;idx++) {
                res[idx] -= t2;
            }

            /* Do multigrid step with solution returned in sg_twovpsi */
            eig_pre[0] = 2;
            eig_post[0] = 2;
            levels=1;
            mgrid_solv (sg_twovpsi, res, work2,
                        dimx, dimy, dimz, hxgrid,
                        hygrid, hzgrid, 0, get_neighbors(), levels, eig_pre, eig_post, 1, 1.0, 0.0,
                        get_NX_GRID(), get_NY_GRID(), get_NZ_GRID(),
                        get_PX_OFFSET(), get_PY_OFFSET(), get_PZ_OFFSET(),
                        get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), ct.boundaryflag);

            for(idx = 0;idx <P0_BASIS;idx++) {
                res[idx] = 0.0;
            }
            pack_stop_axpy (sg_twovpsi, res, 1.0, dimx, dimy, dimz);
            t1 = ct.potential_acceleration_poisson_step;
            if(sp->occupation[0] < 0.5) t1 = 0.0;
#if HYBRID_MODEL
            pthread_mutex_lock(&vtot_sync_mutex);
#endif
            for(idx = 0;idx <P0_BASIS;idx++) {
               vtot_psi[idx] = vtot_psi[idx] + t1 * res[idx];
            }
#if HYBRID_MODEL
            pthread_mutex_unlock(&vtot_sync_mutex);
#endif
        }


    } // end if

    /* Release our memory */
    my_free (nvtot_psi);
    my_free (saved_psi);
    my_free (res2);
#if !BATCH_NLS
    my_free (nv);
    my_free (ns);
#endif
    my_free (work2);
    my_free (sg_twovpsi);
    my_free (work1);
    my_free (res);
    my_free (sg_psi);


}                               /* end mg_eig_state */

#else


/* Complex version */
void mg_eig_state (STATE * sp, int tid, rmg_double_t * vtot_psi)
{

    int idx, cycles;
    int nits, pbasis, sbasis;
    rmg_double_t eig = 0.0, eigR, eigI, diag, t2, t2R, t2I;
    rmg_double_t t1, *kdr, *kdi, *work1R, *work2R, *work1I, *work2I, *gx, *gy, *gz;
    rmg_double_t *tmp_psiR, *nvR, *resR, *sg_psiR, *nsR, *sg_twovpsiR;
    rmg_double_t *tmp_psiI, *nvI, *resI, *sg_psiI, *nsI, *sg_twovpsiI;
    int eig_pre[6] = { 0, 2, 2, 2, 2, 2 };
    int eig_post[6] = { 0, 2, 2, 2, 2, 2 };
    int ione = 1;
    int dimx, dimy, dimz, levels;
    rmg_double_t hxgrid, hygrid, hzgrid, sb_step;
    rmg_double_t time1;


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
    my_malloc (sg_psiR, 10 * (sbasis), rmg_double_t);
    resR = sg_psiR + sbasis;
    work2R = resR + sbasis;
    sg_twovpsiR = work2R + 4 * sbasis;
    work1R = sg_twovpsiR + sbasis;
    nvR = work1R + sbasis;
    nsR = nvR + sbasis;

    my_malloc (sg_psiI, 15 * (sbasis), rmg_double_t);
    resI = sg_psiI + sbasis;
    work2I = resI + sbasis;
    sg_twovpsiI = work2I + 4 * sbasis;
    work1I = sg_twovpsiI + sbasis;
    gx = work1I + sbasis;
    gy = gx + sbasis;
    gz = gy + sbasis;
    kdr = gz + sbasis;
    kdi = kdr + sbasis;
    nvI = kdi + sbasis;
    nsI = nvI + sbasis;


    tmp_psiR = sp->psiR;
    tmp_psiI = sp->psiI;


#if MD_TIMERS
    time1 = my_crtc ();
#endif
    /* Get the non-local operator and S acting on psi (nv and ns, respectively) */
    app_nls (tmp_psiR, tmp_psiI, nvR, nvI, nsR, nsI, pct.oldsintR_local, pct.oldsintI_local, sp->istate, sp->kidx);


#if MD_TIMERS
    rmg_timings (MG_EIG_NLS_TIME, (my_crtc () - time1));
#endif


    /* Smoothing cycles */
    for (cycles = 0; cycles <= nits; cycles++)
    {

#if MD_TIMERS
        time1 = my_crtc ();
#endif
        /* Pack psi into smoothing arrays */
        /*pack_ptos (sg_psiR, tmp_psiR, dimx, dimy, dimz);
        pack_ptos (sg_psiI, tmp_psiI, dimx, dimy, dimz);*/
#if MD_TIMERS
        rmg_timings (MG_EIG_PACK_TIME, (my_crtc () - time1));
#endif


#if MD_TIMERS
        time1 = my_crtc ();
#endif

        /* Apply Mehrstellen left hand operators */
        diag = app_cil_driver (tmp_psiR, work2R, dimx, dimy, dimz,
                              hxgrid, hygrid, hzgrid, ct.kohn_sham_fd_order);


        /* Apply Mehrstellen left and right hand operators */
        diag = app_cil_driver (tmp_psiI, work2I, dimx, dimy, dimz,
                              hxgrid, hygrid, hzgrid, ct.kohn_sham_fd_order);

#if MD_TIMERS
        rmg_timings (MG_EIG_APPCIL_TIME, (my_crtc () - time1));
#endif

        diag = -1.0 / diag;

        /* Apply the gradient operator to psi */
        app_grad (tmp_psiI, (P0_GRID *) gx, (P0_GRID *) gy, (P0_GRID *) gz);

        for (idx = 0; idx <P0_BASIS; idx++)
        {

            kdr[idx] = (ct.kp[sp->kidx].kvec[0] * gx[idx] +
                        ct.kp[sp->kidx].kvec[1] * gy[idx] + ct.kp[sp->kidx].kvec[2] * gz[idx]);

        }

        app_grad (tmp_psiR, (P0_GRID *) gx, (P0_GRID *) gy, (P0_GRID *) gz);
        for (idx = 0; idx <P0_BASIS; idx++)
        {

            kdi[idx] = -(ct.kp[sp->kidx].kvec[0] * gx[idx] +
                         ct.kp[sp->kidx].kvec[1] * gy[idx] + ct.kp[sp->kidx].kvec[2] * gz[idx]);

        }

#if MD_TIMERS
        time1 = my_crtc ();
#endif
        /* Apply Mehrstellen right hand operator on S*psi */
        /*pack_ptos (sg_psiR, nsR, dimx, dimy, dimz);
        pack_ptos (sg_psiI, nsI, dimx, dimy, dimz);*/
#if MD_TIMERS
        rmg_timings (MG_EIG_PACK_TIME, (my_crtc () - time1));
#endif

#if MD_TIMERS
        time1 = my_crtc ();
#endif
        app_cir_driver (nsR, resR, dimx, dimy, dimz, ct.kohn_sham_fd_order);

        /* Apply Mehrstellen right hand operator on imaginary part of S*psi */
        app_cir_driver (nsI, resI, dimx, dimy, dimz, ct.kohn_sham_fd_order);

#if MD_TIMERS
        rmg_timings (MG_EIG_APPCIR_TIME, (my_crtc () - time1));
#endif

#if MD_TIMERS
        time1 = my_crtc ();
#endif
        /* Generate 2 * V * psi */
        genvpsi (tmp_psiR, sg_twovpsiR, vtot_psi, nvR, kdr, ct.kp[sp->kidx].kmag, dimx, dimy, dimz);
        genvpsi (tmp_psiI, sg_twovpsiI, vtot_psi, nvI, kdi, ct.kp[sp->kidx].kmag, dimx, dimy, dimz);
#if MD_TIMERS
        rmg_timings (MG_EIG_GENVPSI_TIME, (my_crtc () - time1));
#endif


#if MD_TIMERS
        time1 = my_crtc ();
#endif

        /* B operating on 2*V*psi stored in work1 */
        app_cir_driver (sg_twovpsiR, work1R, dimx, dimy, dimz, ct.kohn_sham_fd_order);
        app_cir_driver (sg_twovpsiI, work1I, dimx, dimy, dimz, ct.kohn_sham_fd_order);

#if MD_TIMERS
        rmg_timings (MG_EIG_APPCIR_TIME, (my_crtc () - time1));
#endif



        t1 = -ONE;
        QMD_daxpy (pbasis, t1, work2R, ione, work1R, ione);
        QMD_daxpy (pbasis, t1, work2I, ione, work1I, ione);



#if MD_TIMERS
        time1 = my_crtc ();
#endif
        /* If this is the first time through compute the eigenvalue */
        if (cycles == 0)
        {

            eigR = t2R = ZERO;
            for (idx = 0; idx < pbasis; idx++)
            {

                t2R += tmp_psiR[idx] * resR[idx];
                eigR += tmp_psiR[idx] * work1R[idx];

            }

            eigI = t2I = ZERO;
            for (idx = 0; idx < pbasis; idx++)
            {

                t2I += tmp_psiI[idx] * resI[idx];
                eigI += tmp_psiI[idx] * work1I[idx];

            }


            eigR = real_sum_all (eigR, pct.grid_comm);
            t2R = real_sum_all (t2R, pct.grid_comm);
            eigI = real_sum_all (eigI, pct.grid_comm);
            t2I = real_sum_all (t2I, pct.grid_comm);
            sp->eig[0] = (eigR + eigI) / (TWO * (t2R + t2I));
            eig = sp->eig[0];

        }

#if MD_TIMERS
        rmg_timings (MG_EIG_EIGVALUE_TIME, (my_crtc () - time1));
#endif

        /* Next we have to generate the residual vector for smoothing */
        /* or multigridding.                                          */
        t1 = TWO * eig;
        for (idx = 0; idx <P0_BASIS; idx++)
        {

            resR[idx] = t1 * resR[idx] - work1R[idx];
            resI[idx] = t1 * resI[idx] - work1I[idx];

        }

        if (cycles == 0)
        {

            t1 = snrm2 (&pbasis, resR, &ione);
            t1 = t1 * t1;
            t1 = real_sum_all (t1, pct.grid_comm);
            t2 = (rmg_double_t) (ct.psi_nbasis);
            sp->res = ct.hmaxgrid * ct.hmaxgrid * sqrt (t1 / t2) * 0.25;

        }


        /* Now either smooth the wavefunction or do a multigrid cycle */
        if (cycles == ct.eig_parm.gl_pre)
        {

            /* Real components */
            /* Pack the residual data into multigrid array */
#if MD_TIMERS
            time1 = my_crtc ();
#endif
            pack_ptos (sg_psiR, resR, dimx, dimy, dimz);
#if MD_TIMERS
            rmg_timings (MG_EIG_PACK_TIME, (my_crtc () - time1));
#endif

            time1 = my_crtc ();
            trade_images (sg_psiR, dimx, dimy, dimz, get_neighbors(), FULL_FD);
            rmg_timings (MG_EIG_TRADE_TIME, (my_crtc () - time1));

            /* Smooth it once and store the smoothed residual in work1 */
#if MD_TIMERS
            time1 = my_crtc ();
#endif
            app_smooth (sg_psiR, work1R, PX0_GRID, PY0_GRID, PZ0_GRID);

#if MD_TIMERS
            rmg_timings (MG_EIG_APPSMOOTH_TIME, (my_crtc () - time1));
#endif


#if MD_TIMERS
            time1 = my_crtc ();
#endif
            /* Do multigrid step with solution returned in sg_twovpsi */
            mgrid_solv (sg_twovpsiR, work1R, work2R,
                        dimx, dimy, dimz, hxgrid,
                        hygrid, hzgrid, 0, get_neighbors(), levels, eig_pre, eig_post, 1, sb_step);
#if MD_TIMERS
            rmg_timings (MG_EIG_MGRIDSOLV_TIME, (my_crtc () - time1));
#endif


            /* The correction is in a smoothing grid so we use this
             * routine to update the orbital which is stored in a physical grid.
             */
            t1 = -1.0;
            pack_stop_axpy (sg_twovpsiR, tmp_psiR, t1, dimx, dimy, dimz);


            /* Imaginary components */
            /* Pack the residual data into multigrid array */
#if MD_TIMERS
            time1 = my_crtc ();
#endif
            pack_ptos (sg_psiI, resI, dimx, dimy, dimz);
#if MD_TIMERS
            rmg_timings (MG_EIG_PACK_TIME, (my_crtc () - time1));
#endif

#if MD_TIMERS
            time1 = my_crtc ();
#endif
            trade_images (sg_psiI, dimx, dimy, dimz, get_neighbors(), FULL_FD);
#if MD_TIMERS
            rmg_timings (MG_EIG_TRADE_TIME, (my_crtc () - time1));
#endif

            /* Smooth it once and store the smoothed residual in work1 */
            t1 = 145.0;

#if MD_TIMERS
            time1 = my_crtc ();
#endif
            app_smooth (sg_psiI, work1I, PX0_GRID, PY0_GRID, PZ0_GRID);

#if MD_TIMERS
            rmg_timings (MG_EIG_APPSMOOTH_TIME, (my_crtc () - time1));
#endif


#if MD_TIMERS
            time1 = my_crtc ();
#endif
            /* Do multigrid step with solution returned in sg_twovpsi */
            mgrid_solv (sg_twovpsiI, work1I, work2I,
                        dimx, dimy, dimz, hxgrid,
                        hygrid, hzgrid, 0, get_neighbors(), levels, eig_pre, eig_post, 1, sb_step);
#if MD_TIMERS
            rmg_timings (MG_EIG_MGRIDSOLV_TIME, (my_crtc () - time1));
#endif


            /* The correction is in a smoothing grid so we use this
             * routine to update the orbital which is stored in a physical grid.
             */
            t1 = -1.0;
#if MD_TIMERS
            time1 = my_crtc ();
#endif
            pack_stop_axpy (sg_twovpsiI, tmp_psiI, t1, dimx, dimy, dimz);
#if MD_TIMERS
            rmg_timings (MG_EIG_PACK_TIME, (my_crtc () - time1));
#endif



        }
        else
        {


            t1 = ct.eig_parm.gl_step * diag;
/*            if(cycles == 0)t1 = 2.0 * diag;
            if(cycles == 1)t1 = 1.0 * diag;*/

            /* Update wavefuntion */
            QMD_daxpy (pbasis, t1, resR, ione, tmp_psiR, ione);
            QMD_daxpy (pbasis, t1, resI, ione, tmp_psiI, ione);


        }

    }                           /* end for */



    /* Release our memory */
    my_free (sg_psiI);
    my_free (sg_psiR);


}                               /* end mg_eig_state */


#endif




/******/

