
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "common_prototypes.h"
#include "BaseThread.h"

#if GAMMA_PT
#include "hybrid.h"
void subdiag_app_A_one (STATE *sp, rmg_double_t * a_psi, rmg_double_t * s_psi, rmg_double_t * vtot_eig);
void subdiag_app_B_one (STATE *sp, rmg_double_t * b_psi);


/*Applies A operator to all wavefunctions*/
void subdiag_app_A (STATE * states, rmg_double_t * a_psi, rmg_double_t * s_psi, rmg_double_t * vtot_eig)
{
    int istate, st1, ist, istop, P0_BASIS;
    STATE *sp;
    rmg_double_t time1, time2;

#if BATCH_NLS
    time1 = my_crtc();
    app_nls_batch (states, pct.nv, s_psi, pct.Bns, pct.newsintR_local);
#endif

    P0_BASIS = get_P0_BASIS();

    // Each thread applies the operator to one wavefunction
    istop = ct.num_states / ct.THREADS_PER_NODE;
    istop = istop * ct.THREADS_PER_NODE;     
    for(st1=0;st1 < istop;st1+=ct.THREADS_PER_NODE) {
        SCF_THREAD_CONTROL thread_control[MAX_RMG_THREADS];
        for(ist = 0;ist < ct.THREADS_PER_NODE;ist++) {
            thread_control[ist].job = HYBRID_SUBDIAG_APP_A;
            thread_control[ist].sp = &states[st1 + ist];
            thread_control[ist].p1 = &a_psi[(st1 + ist) * P0_BASIS];
            thread_control[ist].p2 = &s_psi[(st1 + ist) * P0_BASIS];
            thread_control[ist].vtot = vtot_eig;
            set_pptr(ist, &thread_control[ist]);
        }

        // Thread tasks are set up so wake them
        run_thread_tasks(ct.THREADS_PER_NODE);

    }

    // Process any remaining orbitals serially
    for(st1 = istop;st1 < ct.num_states;st1++) {
        subdiag_app_A_one (&states[st1], &a_psi[st1 * P0_BASIS], &s_psi[st1 *  P0_BASIS], vtot_eig);
    }
#if GPU_ENABLED
    cuCtxSynchronize();
#endif
}

// Applies A operator to one wavefunction
void subdiag_app_A_one (STATE *sp, rmg_double_t * a_psi, rmg_double_t * s_psi, rmg_double_t * vtot_eig)
{
    int kidx, idx, istate, sbasis, tid, P0_BASIS;
    rmg_double_t *sg_twovpsi, *tmp_psi, *work2, *work1, *work3;
#    if MD_TIMERS
    rmg_double_t time1;
#    endif

    P0_BASIS = get_P0_BASIS();

    tid = get_thread_tid();
    if(tid < 0) tid = 0;  // OK in this case

    sbasis = sp->sbasis;

    my_malloc (work3, sbasis, rmg_double_t);

#if !BATCH_NLS
    my_malloc (work2, sbasis, rmg_double_t);
#endif
    my_malloc (sg_twovpsi, sbasis, rmg_double_t);
    kidx = 0;

    work1 = a_psi;

    tmp_psi = sp->psiR;

#   if MD_TIMERS
    time1 = my_crtc ();
#   endif
    /* A operating on psi stored in work3 */
    app_cil_driver (tmp_psi, work3, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_hxgrid(),
                   get_hygrid(), get_hzgrid(), ct.kohn_sham_fd_order);

#   if MD_TIMERS
#   endif


#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* Apply non-local operator to psi and store in work2 */
#if !BATCH_NLS
    app_nls (tmp_psi, NULL, work2, NULL, s_psi, NULL, pct.newsintR_local, NULL, sp->istate, kidx);
#   if MD_TIMERS
#   endif

#else
    work2 = &pct.nv[sp->istate *  P0_BASIS];
#endif



#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* Generate 2*V*psi and store it in a smoothing grid and store in sg_twovpsi */
        if((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0)) {
            if(ct.scf_steps == 0) {
                    genvpsi (tmp_psi, sg_twovpsi, vtot_eig, work2, 0.0, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID());
            }
            else {
                    genvpsi (tmp_psi, sg_twovpsi, sp->dvhxc, work2, 0.0, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID());
            }
        }
        else {
            genvpsi (tmp_psi, sg_twovpsi, vtot_eig, work2, 0.0, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID());
        }


#   if MD_TIMERS
#   endif

#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* B operating on 2*V*psi stored in work1 */
    app_cir_driver (sg_twovpsi, work1, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), ct.kohn_sham_fd_order);
    for(idx = 0; idx <  P0_BASIS; idx++) work1[idx] += TWO * work2[idx];
#   if MD_TIMERS
#   endif

    /* Pack psi into smoothing array */
    //pack_ptos (sg_psi, tmp_psi, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID());




    for (idx = 0; idx < P0_BASIS; idx++)
        work1[idx] = 0.5 * get_vel() * (work1[idx] - work3[idx]);

#if !BATCH_NLS
    my_free (work2);
#endif
    my_free(sg_twovpsi);
    my_free(work3);

}                               /* subdiag_app_A_one */


void subdiag_app_B (STATE * states, rmg_double_t * b_psi)
{
    int st1, ist, istate, istop, P0_BASIS;
    STATE *sp;

    P0_BASIS = get_P0_BASIS();

    // Each thread applies the operator to one wavefunction
    istop = ct.num_states / ct.THREADS_PER_NODE;
    istop = istop * ct.THREADS_PER_NODE;
    for(st1=0;st1 < istop;st1+=ct.THREADS_PER_NODE) {
        SCF_THREAD_CONTROL thread_control[MAX_RMG_THREADS];
        for(ist = 0;ist < ct.THREADS_PER_NODE;ist++) {
            thread_control[ist].job = HYBRID_SUBDIAG_APP_B;
            thread_control[ist].sp = &states[st1 + ist];
            thread_control[ist].p1 = &b_psi[(st1 + ist) * P0_BASIS];
            set_pptr(ist, &thread_control[ist]);
        }

        // Thread tasks are set up so wake them
        run_thread_tasks(ct.THREADS_PER_NODE);

    }

    // Process any remaining orbitals serially
    for(st1 = istop;st1 < ct.num_states;st1++) {
        subdiag_app_B_one(&states[st1], &b_psi[st1 * P0_BASIS]);
    }

#if GPU_ENABLED
    cuCtxSynchronize();
#endif

}


void subdiag_app_B_one (STATE *sp, rmg_double_t * b_psi)
{
    int istate, pbasis, ione=1;
    rmg_double_t *work2, *work1;

#    if MD_TIMERS
    rmg_double_t time1;
#    endif

    pbasis = sp->pbasis;

    my_malloc (work2, pbasis, rmg_double_t);

    work1 = b_psi;

    /*Pack S|psi> into smoothing array */
    //pack_ptos (sg_psi, work1, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID());
    scopy (&pbasis, work1, &ione, work2, &ione);


#   if MD_TIMERS
        time1 = my_crtc ();
#   endif
    /*B operating on S|psi> and store in work3 */
    app_cir_driver (work2, work1, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), ct.kohn_sham_fd_order);
#   if MD_TIMERS
#   endif

    my_free (work2);

}                               /* subdiag_app_B_one */

#endif
