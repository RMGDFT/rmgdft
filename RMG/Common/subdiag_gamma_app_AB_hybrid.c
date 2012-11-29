
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

    enter_threaded_region();
    scf_barrier_init(THREADS_PER_NODE);

    // Each thread applies the operator to one wavefunction
    istop = ct.num_states / THREADS_PER_NODE;
    istop = istop * THREADS_PER_NODE;     
    for(st1=0;st1 < istop;st1+=THREADS_PER_NODE) {
        for(ist = 0;ist < THREADS_PER_NODE;ist++) {
            thread_control[ist].job = HYBRID_SUBDIAG_APP_A;
            thread_control[ist].sp = &states[st1 + ist];
            thread_control[ist].p1 = &a_psi[(st1 + ist) * P0_BASIS];
            thread_control[ist].p2 = &s_psi[(st1 + ist) * P0_BASIS];
            thread_control[ist].vtot = vtot_eig;
        }

        // Thread tasks are set up so wake them
        wake_threads(THREADS_PER_NODE);

        // Then wait for them to finish this task
        wait_for_threads(THREADS_PER_NODE);

    }

    scf_barrier_destroy();
    leave_threaded_region();

    // Process any remaining orbitals serially
    for(st1 = istop;st1 < ct.num_states;st1++) {
        subdiag_app_A_one (&states[st1], &a_psi[st1 * P0_BASIS], &s_psi[st1 * P0_BASIS], vtot_eig);
    }
}

// Applies A operator to one wavefunction
void subdiag_app_A_one (STATE *sp, REAL * a_psi, REAL * s_psi, REAL * vtot_eig)
{
    int kidx, idx, istate, sbasis;
    REAL *sg_twovpsi, *tmp_psi, *work2, *work1;
#    if MD_TIMERS
    REAL time1;
#    endif


    sbasis = sp->sbasis;
    my_malloc (work2, 2 * sbasis, REAL);
    sg_twovpsi = work2 + sbasis;
    kidx = 0;

    work1 = a_psi;

    tmp_psi = sp->psiR;


#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* Apply non-local operator to psi and store in work2 */
    app_nls (tmp_psi, NULL, work2, NULL, s_psi, NULL, pct.newsintR_local, NULL, sp->istate, kidx);
#   if MD_TIMERS
    rmg_timings (DIAG_NL_TIME, (my_crtc () - time1));
#   endif



#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* Generate 2*V*psi and store it in a smoothing grid and store in sg_twovpsi */
        if(ct.potential_acceleration) {
            if(ct.scf_steps == 0) {
                    genvpsi (tmp_psi, sg_twovpsi, vtot_eig, work2, NULL, 0.0, PX0_GRID, PY0_GRID, PZ0_GRID);
            }
            else {
                    genvpsi (tmp_psi, sg_twovpsi, sp->dvhxc, work2, NULL, 0.0, PX0_GRID, PY0_GRID, PZ0_GRID);
            }
        }
        else {
            genvpsi (tmp_psi, sg_twovpsi, vtot_eig, work2, NULL, 0.0, PX0_GRID, PY0_GRID, PZ0_GRID);
        }


#   if MD_TIMERS
        rmg_timings (DIAG_GENVPSI_TIME, (my_crtc () - time1));
#   endif

#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* B operating on 2*V*psi stored in work1 */
    app_cir_driver (sg_twovpsi, work1, PX0_GRID, PY0_GRID, PZ0_GRID, ct.kohn_sham_fd_order);
#   if MD_TIMERS
       rmg_timings (DIAG_APPCIR_TIME, (my_crtc () - time1));
#   endif

    /* Pack psi into smoothing array */
    //pack_ptos (sg_psi, tmp_psi, PX0_GRID, PY0_GRID, PZ0_GRID);


#   if MD_TIMERS
        time1 = my_crtc ();
#   endif
    /* A operating on psi stored in work2 */
    app_cil_driver (tmp_psi, work2, PX0_GRID, PY0_GRID, PZ0_GRID, sp->hxgrid,
                   sp->hygrid, sp->hzgrid, ct.kohn_sham_fd_order);

#   if MD_TIMERS
        rmg_timings (DIAG_APPCIL_TIME, (my_crtc () - time1));
#   endif

    for (idx = 0; idx < P0_BASIS; idx++)
        work1[idx] = 0.5 * ct.vel * (work1[idx] - work2[idx]);

    my_free (work2);

}                               /* subdiag_app_A_one */


void subdiag_app_B (STATE * states, REAL * b_psi)
{
    int st1, ist, istate, istop;
    STATE *sp;

    enter_threaded_region();
    scf_barrier_init(THREADS_PER_NODE);

    // Each thread applies the operator to one wavefunction
    istop = ct.num_states / THREADS_PER_NODE;
    istop = istop * THREADS_PER_NODE;
    for(st1=0;st1 < istop;st1+=THREADS_PER_NODE) {
        for(ist = 0;ist < THREADS_PER_NODE;ist++) {
            thread_control[ist].job = HYBRID_SUBDIAG_APP_B;
            thread_control[ist].sp = &states[st1 + ist];
            thread_control[ist].p1 = &b_psi[(st1 + ist) * P0_BASIS];
        }

        // Thread tasks are set up so wake them
        wake_threads(THREADS_PER_NODE);

        // Then wait for them to finish this task
        wait_for_threads(THREADS_PER_NODE);

    }

    scf_barrier_destroy();
    leave_threaded_region();

    // Process any remaining orbitals serially
    for(st1 = istop;st1 < ct.num_states;st1++) {
        subdiag_app_B_one(&states[st1], &b_psi[st1 * P0_BASIS]);
    }

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
    //pack_ptos (sg_psi, work1, PX0_GRID, PY0_GRID, PZ0_GRID);
    scopy (&pbasis, work1, &ione, work2, &ione);


#   if MD_TIMERS
        time1 = my_crtc ();
#   endif
    /*B operating on S|psi> and store in work3 */
    app_cir_driver (work2, work1, PX0_GRID, PY0_GRID, PZ0_GRID, ct.kohn_sham_fd_order);
#   if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME2, (my_crtc () - time1));
#   endif

    my_free (work2);

}                               /* subdiag_app_B_one */

#endif
#endif
