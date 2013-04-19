
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

#if HYBRID_MODEL
#if GAMMA_PT
#include "hybrid.h"
#include <pthread.h>
void subdiag_app_A_one (STATE *sp, REAL * a_psi, REAL * s_psi, REAL * vtot_eig);
void subdiag_app_B_one (STATE *sp, REAL * b_psi);


/*Applies A operator to all wavefunctions*/
void subdiag_app_A (STATE * states, REAL * a_psi, REAL * s_psi, REAL * vtot_eig)
{
    int istate, st1, ist, istop;
    STATE *sp;

#if BATCH_NLS
    app_nls_batch (states, pct.nv, s_psi, pct.newsintR_local);
#endif

    enter_threaded_region();
    scf_barrier_init(ct.THREADS_PER_NODE);

    // Each thread applies the operator to one wavefunction
    istop = ct.num_states / ct.THREADS_PER_NODE;
    istop = istop * ct.THREADS_PER_NODE;     
    for(st1=0;st1 < istop;st1+=ct.THREADS_PER_NODE) {
        for(ist = 0;ist < ct.THREADS_PER_NODE;ist++) {
            thread_control[ist].job = HYBRID_SUBDIAG_APP_A;
            thread_control[ist].sp = &states[st1 + ist];
            thread_control[ist].p1 = &a_psi[(st1 + ist) *pct.P0_BASIS];
            thread_control[ist].p2 = &s_psi[(st1 + ist) *pct.P0_BASIS];
            thread_control[ist].vtot = vtot_eig;
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
        subdiag_app_A_one (&states[st1], &a_psi[st1 *pct.P0_BASIS], &s_psi[st1 * pct.P0_BASIS], vtot_eig);
    }
#if GPU_ENABLED
    cuCtxSynchronize();
#endif
}

// Applies A operator to one wavefunction
void subdiag_app_A_one (STATE *sp, REAL * a_psi, REAL * s_psi, REAL * vtot_eig)
{
    int kidx, idx, istate, sbasis, tid;
    REAL *sg_twovpsi, *tmp_psi, *work2, *work1, *work3;
#    if MD_TIMERS
    REAL time1;
#    endif

#if HYBRID_MODEL
    tid = get_thread_tid();
    if(tid < 0) tid = 0;  // OK in this case
#else
    tid = 0;
#endif

    sbasis = sp->sbasis;

#if GPU_FD_ENABLED
    cudaStream_t *cstream;
    cstream = get_thread_cstream();
    work3 = &ct.gpu_host_temp3[tid * sbasis];
#else
    my_malloc (work3, sbasis, REAL);
#endif

#if !BATCH_NLS
    my_malloc (work2, sbasis, REAL);
#endif
    my_malloc (sg_twovpsi, sbasis, REAL);
    kidx = 0;

    work1 = a_psi;

    tmp_psi = sp->psiR;

#   if MD_TIMERS
    time1 = my_crtc ();
#   endif
    /* A operating on psi stored in work3 */
    app_cil_driver (tmp_psi, work3, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, sp->hxgrid,
                   sp->hygrid, sp->hzgrid, ct.kohn_sham_fd_order);

#   if MD_TIMERS
    rmg_timings (DIAG_APPCIL_TIME, (my_crtc () - time1));
#   endif


#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* Apply non-local operator to psi and store in work2 */
#if !BATCH_NLS
    app_nls (tmp_psi, NULL, work2, NULL, s_psi, NULL, pct.newsintR_local, NULL, sp->istate, kidx);
#else
    work2 = &pct.nv[sp->istate * pct.P0_BASIS];
#endif

#   if MD_TIMERS
    rmg_timings (DIAG_NL_TIME, (my_crtc () - time1));
#   endif



#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* Generate 2*V*psi and store it in a smoothing grid and store in sg_twovpsi */
        if((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0)) {
            if(ct.scf_steps == 0) {
                    genvpsi (tmp_psi, sg_twovpsi, vtot_eig, work2, NULL, 0.0, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);
            }
            else {
                    genvpsi (tmp_psi, sg_twovpsi, sp->dvhxc, work2, NULL, 0.0, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);
            }
        }
        else {
            genvpsi (tmp_psi, sg_twovpsi, vtot_eig, work2, NULL, 0.0, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);
        }


#   if MD_TIMERS
        rmg_timings (DIAG_GENVPSI_TIME, (my_crtc () - time1));
#   endif

#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* B operating on 2*V*psi stored in work1 */
    app_cir_driver (sg_twovpsi, work1, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, ct.kohn_sham_fd_order);
#   if MD_TIMERS
       rmg_timings (DIAG_APPCIR_TIME, (my_crtc () - time1));
#   endif

    /* Pack psi into smoothing array */
    //pack_ptos (sg_psi, tmp_psi, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);




#if GPU_FD_ENABLED
    cuStreamSynchronize(*cstream);
#endif
    for (idx = 0; idx <pct.P0_BASIS; idx++)
        work1[idx] = 0.5 * ct.vel * (work1[idx] - work3[idx]);

#if !BATCH_NLS
    my_free (work2);
#endif
    my_free(sg_twovpsi);

#if !GPU_FD_ENABLED
    my_free(work3);
#endif

}                               /* subdiag_app_A_one */


void subdiag_app_B (STATE * states, REAL * b_psi)
{
    int st1, ist, istate, istop;
    STATE *sp;

    enter_threaded_region();
    scf_barrier_init(ct.THREADS_PER_NODE);

    // Each thread applies the operator to one wavefunction
    istop = ct.num_states / ct.THREADS_PER_NODE;
    istop = istop * ct.THREADS_PER_NODE;
    for(st1=0;st1 < istop;st1+=ct.THREADS_PER_NODE) {
        for(ist = 0;ist < ct.THREADS_PER_NODE;ist++) {
            thread_control[ist].job = HYBRID_SUBDIAG_APP_B;
            thread_control[ist].sp = &states[st1 + ist];
            thread_control[ist].p1 = &b_psi[(st1 + ist) *pct.P0_BASIS];
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
        subdiag_app_B_one(&states[st1], &b_psi[st1 *pct.P0_BASIS]);
    }

#if GPU_ENABLED
    cuCtxSynchronize();
#endif

}


void subdiag_app_B_one (STATE *sp, REAL * b_psi)
{
    int istate, pbasis, ione=1;
    REAL *work2, *work1;

#    if MD_TIMERS
    REAL time1;
#    endif

    pbasis = sp->pbasis;

    my_malloc (work2, pbasis, REAL);

    work1 = b_psi;

    /*Pack S|psi> into smoothing array */
    //pack_ptos (sg_psi, work1, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);
    scopy (&pbasis, work1, &ione, work2, &ione);


#   if MD_TIMERS
        time1 = my_crtc ();
#   endif
    /*B operating on S|psi> and store in work3 */
    app_cir_driver (work2, work1, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, ct.kohn_sham_fd_order);
#   if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME2, (my_crtc () - time1));
#   endif

    my_free (work2);

}                               /* subdiag_app_B_one */

#endif
#endif
