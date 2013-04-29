
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

#if HYBRID_MODEL
#if GAMMA_PT
#include "hybrid.h"
#include <pthread.h>
void subdiag_app_AB_one (STATE *sp, REAL * a_psi, REAL * b_psi, REAL * vtot_eig);


/*Applies A operator to all wavefunctions*/
void subdiag_app_AB (STATE * states, REAL * a_psi, REAL * b_psi, REAL * vtot_eig)
{
    int istate, st1, ist, istop;
    STATE *sp;
    REAL time1, time2;
    REAL *vtot_eig_s;
    int ione = 1;
#if !BATCH_NLS
    error_handler("set BATCJ_NLS 1, other mode not programmed here");
#endif
    
    int idx, dimx = pct.PX0_GRID, dimy = pct.PY0_GRID, dimz = pct.PZ0_GRID;
    idx = (pct.PX0_GRID + 4) * (pct.PY0_GRID + 4) * (pct.PZ0_GRID + 4) ;
    my_malloc(vtot_eig_s, idx, REAL);

// trade images for vtot_eig and stored in vtot_eig_s. It will be used
// in app_cilr
//
    if(ct.kohn_sham_fd_order == APP_CI_FOURTH) 
        trade_imagesx (vtot_eig, vtot_eig_s, dimx, dimy, dimz, 1, FULL_FD);

    if(ct.kohn_sham_fd_order == APP_CI_SIXTH) 
        trade_imagesx (vtot_eig, vtot_eig_s, dimx, dimy, dimz, 2, FULL_FD);

    time1 = my_crtc();
    app_nls_batch (states, pct.nv, pct.ns, pct.Bns, pct.newsintR_local);
    rmg_timings (DIAG_NL_TIME, (my_crtc () - time1));

    enter_threaded_region();
    scf_barrier_init(ct.THREADS_PER_NODE);

    // Each thread applies the operator to one wavefunction
    istop = ct.num_states / ct.THREADS_PER_NODE;
    istop = istop * ct.THREADS_PER_NODE;     
    for(st1=0;st1 < istop;st1+=ct.THREADS_PER_NODE) {
        for(ist = 0;ist < ct.THREADS_PER_NODE;ist++) {
            thread_control[ist].job = HYBRID_SUBDIAG_APP_AB;
            thread_control[ist].sp = &states[st1 + ist];
            thread_control[ist].p1 = &a_psi[(st1 + ist) *pct.P0_BASIS];
            thread_control[ist].p2 = &b_psi[(st1 + ist) *pct.P0_BASIS];
            thread_control[ist].vtot = vtot_eig_s;
        }

        // Thread tasks are set up so wake them
        wake_threads(ct.THREADS_PER_NODE);

        // Then wait for them to finish this task
        wait_for_threads(ct.THREADS_PER_NODE);

    }

    scf_barrier_destroy();
    leave_threaded_region();

    // Process any remaining orbitals serially
    for(st1 = istop;st1 < ct.num_states;st1++) {
        subdiag_app_AB_one (&states[st1], &a_psi[st1 *pct.P0_BASIS], &b_psi[st1 * pct.P0_BASIS], vtot_eig_s);
    }
#if GPU_ENABLED
    cuCtxSynchronize();
#endif

    my_free(vtot_eig_s);
}

// Applies A operator to one wavefunction
void subdiag_app_AB_one (STATE *sp, REAL * a_psi, REAL * b_psi, REAL * vtot_eig_s)
{
    int idx, istate, sbasis, tid;
    REAL *sg_twovpsi, *tmp_psi, *work2;
    int dimx = pct.PX0_GRID, dimy = pct.PY0_GRID, dimz = pct.PZ0_GRID;
    int ione = 1;
#    if MD_TIMERS
    REAL time1;
#    endif


    tmp_psi = sp->psiR;

#   if MD_TIMERS
    time1 = my_crtc ();
#   endif
    /* Generate 2*V*psi and store it in a smoothing grid and store in sg_twovpsi */
    if((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0)) {
        if(ct.scf_steps > 0) 
        {
            if(ct.kohn_sham_fd_order == APP_CI_FOURTH) 
                trade_imagesx (sp->dvhxc, vtot_eig_s, dimx, dimy, dimz, 1, FULL_FD);

            if(ct.kohn_sham_fd_order == APP_CI_SIXTH) 
                trade_imagesx (sp->dvhxc, vtot_eig_s, dimx, dimy, dimz, 2, FULL_FD);
        }
    }

    /* A operating on psi stored in work3 */
    app_cilr_driver (tmp_psi, a_psi, b_psi, vtot_eig_s, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, sp->hxgrid,
            sp->hygrid, sp->hzgrid, ct.kohn_sham_fd_order);

#   if MD_TIMERS
    rmg_timings (DIAG_APPCIL_TIME, (my_crtc () - time1));
#   endif


    /* Apply non-local operator to psi and store in work2 */
    work2 = &pct.nv[sp->istate * pct.P0_BASIS];


    for(idx = 0; idx < pct.P0_BASIS; idx++) a_psi[idx] += TWO * work2[idx];

    for (idx = 0; idx <pct.P0_BASIS; idx++)
        a_psi[idx] = 0.5 * ct.vel * a_psi[idx];

    work2 = &pct.Bns[sp->istate * pct.P0_BASIS];
    for(idx = 0; idx < pct.P0_BASIS; idx++) b_psi[idx] += work2[idx];


}                               /* subdiag_app_A_one */


#endif
#endif
