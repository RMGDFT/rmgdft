

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

#if !HYBRID_MODEL


#  if GAMMA_PT
/*Applies A operator to all wavefunctions*/
void subdiag_app_A (STATE * states, rmg_double_t * a_psi, rmg_double_t * s_psi, rmg_double_t * vtot_eig)
{
    int i, kidx, idx, istate, sbasis;
    rmg_double_t *sg_twovpsi, *tmp_psi, *work2, *work1;
    STATE *sp;
#    if MD_TIMERS
    rmg_double_t time1;
#    endif

#if BATCH_NLS
    time1 = my_crtc();
    app_nls_batch (states, pct.nv, s_psi, pct.Bns, pct.newsintR_local);
    rmg_timings (DIAG_NL_TIME, (my_crtc () - time1));
#endif

    sbasis = states[0].sbasis;
#if !BATCH_NLS
    my_malloc (work2, sbasis, rmg_double_t);
#endif
    my_malloc (sg_twovpsi, sbasis, rmg_double_t);
    kidx = 0;

    work1 = a_psi;


    for (istate = 0; istate < ct.num_states; istate++)
    {

        sp = &states[istate];
        tmp_psi = sp->psiR;


#if !BATCH_NLS
    app_nls (tmp_psi, NULL, work2, NULL, s_psi, NULL, pct.newsintR_local, NULL, sp->istate, kidx);
#   if MD_TIMERS
    rmg_timings (DIAG_NL_TIME, (my_crtc () - time1));
#   endif

#else
    work2 = &pct.nv[sp->istate * get_P0_BASIS()];
#endif


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif

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
    

#    if MD_TIMERS
        rmg_timings (DIAG_GENVPSI_TIME, (my_crtc () - time1));
#    endif


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* B operating on 2*V*psi stored in work1 */
        app_cir_driver (sg_twovpsi, work1, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, ct.kohn_sham_fd_order);
    for(i = 0; i < get_P0_BASIS(); i++) work1[i] += work2[i];
#    if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME, (my_crtc () - time1));
#    endif


        /* Pack psi into smoothing array */
        //pack_ptos (sg_psi, tmp_psi, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* A operating on psi stored in work2 */
        app_cil_driver (tmp_psi, work2, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, sp->hxgrid,
                       sp->hygrid, sp->hzgrid, ct.kohn_sham_fd_order);

#    if MD_TIMERS
        rmg_timings (DIAG_APPCIL_TIME, (my_crtc () - time1));
#    endif

        for (idx = 0; idx < get_P0_BASIS(); idx++)
            work1[idx] = 0.5 * ct.vel * (work1[idx] - work2[idx]);




        work1 += get_P0_BASIS();
        s_psi += get_P0_BASIS();
    }

    my_free (sg_twovpsi);
#if !BATCH_NLS
    my_free (work2);
#endif

}                               /* subdiag_app_A */


/*Applies B operator to all wavefunctions*/
/*On input b_psi contains s-operator applied to wavefunction*/
void subdiag_app_B (STATE * states, rmg_double_t * b_psi)
{
    int istate, pbasis, ione=1;
    rmg_double_t *work2, *work1;
#    if MD_TIMERS
    rmg_double_t time1;
#    endif

    pbasis = states[0].pbasis;


    my_malloc (work2, pbasis, rmg_double_t);

    work1 = b_psi;


    for (istate = 0; istate < ct.num_states; istate++)
    {

        /*Pack S|psi> into smoothing array */
        //pack_ptos (sg_psi, work1, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);
	scopy (&pbasis, work1, &ione, work2, &ione);


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /*B operating on S|psi> and store in work3 */
        app_cir_driver (work2, work1, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, ct.kohn_sham_fd_order);
#    if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME2, (my_crtc () - time1));
#    endif


        work1 += get_P0_BASIS();
    }

    my_free (work2);

}                               /* subdiag_app_B */


/*Applies A operator to all wavefunctions*/
void subdiag_app_AB (STATE * states, rmg_double_t * a_psi, rmg_double_t * b_psi, rmg_double_t * vtot_eig)
{
    int istate, st1, ist, istop;
    STATE *sp;
    rmg_double_t time1, time2;
    rmg_double_t *vtot_eig_s;
    int ione = 1;
#if !BATCH_NLS
    error_handler("set BATCJ_NLS 1, other mode not programmed here");
#endif

    int idx, dimx = pct.PX0_GRID, dimy = pct.PY0_GRID, dimz = pct.PZ0_GRID;
    idx = (pct.PX0_GRID + 4) * (pct.PY0_GRID + 4) * (pct.PZ0_GRID + 4) ;
    my_malloc(vtot_eig_s, idx, rmg_double_t);

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


    // Process any remaining orbitals serially
    for(st1 = 0;st1 < ct.num_states;st1++) {
        subdiag_app_AB_one (&states[st1], &a_psi[st1 * get_P0_BASIS()], &b_psi[st1 * get_P0_BASIS()], vtot_eig_s);
    }
#if GPU_ENABLED
    cuCtxSynchronize();
#endif

    my_free(vtot_eig_s);
}

// Applies A operator to one wavefunction
void subdiag_app_AB_one (STATE *sp, rmg_double_t * a_psi, rmg_double_t * b_psi, rmg_double_t * vtot_eig_s)
{
    int idx, istate, sbasis, tid;
    rmg_double_t *sg_twovpsi, *tmp_psi, *work2;
    int dimx = pct.PX0_GRID, dimy = pct.PY0_GRID, dimz = pct.PZ0_GRID;
    int ione = 1;
#    if MD_TIMERS
    rmg_double_t time1;
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
    work2 = &pct.nv[sp->istate * get_P0_BASIS()];


    for(idx = 0; idx < get_P0_BASIS(); idx++) a_psi[idx] += TWO * work2[idx];

    for (idx = 0; idx < get_P0_BASIS(); idx++)
        a_psi[idx] = 0.5 * ct.vel * a_psi[idx];

    work2 = &pct.Bns[sp->istate * get_P0_BASIS()];
    for(idx = 0; idx < get_P0_BASIS(); idx++) b_psi[idx] += work2[idx];


}                               /* subdiag_app_A_one */


#else

/*Applies A operator to all wavefunctions*/
void subdiag_app_A (STATE * states, rmg_double_t * a_psiR, rmg_double_t * a_psiI, rmg_double_t * s_psiR, rmg_double_t * s_psiI, rmg_double_t * vtot_eig)
{
    int kidx, idx, istate, sbasis;
    rmg_double_t *sg_twovpsiR, *sg_twovpsiI, *tmp_psiR, *tmp_psiI, *work2R, *work2I,
        *work1R, *work1I;
    rmg_double_t *gx, *gy, *gz, *kdr;
    STATE *sp;
#    if MD_TIMERS
    rmg_double_t time1;
#    endif



    sbasis = states[0].sbasis;
    my_malloc (work2R, 8 * sbasis, rmg_double_t);
    work2I = work2R + sbasis;
    sg_twovpsiR = work2I + sbasis;
    sg_twovpsiI = sg_twovpsiR + sbasis;
    gx = sg_twovpsiI + sbasis;
    gy = gx + sbasis;
    gz = gy + sbasis;
    kdr = gz + sbasis;



    kidx = states[0].kidx;

    work1R = a_psiR;
    work1I = a_psiI;




    for (istate = 0; istate < ct.num_states; istate++)
    {

        sp = &states[istate];
        tmp_psiR = sp->psiR;
        tmp_psiI = sp->psiI;


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* Apply non-local operator to psi and store in work2 */
        app_nls (tmp_psiR, tmp_psiI, work2R, work2I, s_psiR, s_psiI, pct.newsintR_local, pct.newsintI_local, FALSE, kidx);
#    if MD_TIMERS
        rmg_timings (DIAG_NL_TIME, (my_crtc () - time1));
#    endif



        /* Pack psi into smoothing array and get image data */
        /*pack_ptos (sg_psiR, tmp_psiR, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);
        pack_ptos (sg_psiI, tmp_psiI, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);*/


        /* Apply the gradient operator to psi */
        app_grad (tmp_psiI, (P0_GRID *) gx, (P0_GRID *) gy, (P0_GRID *) gz);
        for (idx = 0; idx < get_P0_BASIS(); idx++)
            kdr[idx] = (ct.kp[sp->kidx].kvec[0] * gx[idx] +
                        ct.kp[sp->kidx].kvec[1] * gy[idx] + ct.kp[sp->kidx].kvec[2] * gz[idx]);




#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* Generate 2*V*psi and store in sg_twovpsi */
        genvpsi (tmp_psiR, sg_twovpsiR, vtot_eig, work2R, kdr, ct.kp[sp->kidx].kmag, pct.PX0_GRID,
                 pct.PY0_GRID, pct.PZ0_GRID);
#    if MD_TIMERS
        rmg_timings (DIAG_GENVPSI_TIME, (my_crtc () - time1));
#    endif


        /* Apply the gradient operator to psi */
        app_grad (tmp_psiR, (P0_GRID *) gx, (P0_GRID *) gy, (P0_GRID *) gz);
        for (idx = 0; idx < get_P0_BASIS(); idx++)
            kdr[idx] = -(ct.kp[sp->kidx].kvec[0] * gx[idx] +
                         ct.kp[sp->kidx].kvec[1] * gy[idx] + ct.kp[sp->kidx].kvec[2] * gz[idx]);


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* Generate 2 * V * psiI and store it in a smoothing grid and store in sg_twovpsiI */
        genvpsi (tmp_psiI, sg_twovpsiI, vtot_eig, work2I, kdr, ct.kp[sp->kidx].kmag, pct.PX0_GRID,
                 pct.PY0_GRID, pct.PZ0_GRID);
#    if MD_TIMERS
        rmg_timings (DIAG_GENVPSI_TIME, (my_crtc () - time1));
#    endif




#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* B operating on 2*V*psi stored in work1 */
        app_cir_driver (sg_twovpsiR, work1R, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, ct.kohn_sham_fd_order);
        app_cir_driver (sg_twovpsiI, work1I, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, ct.kohn_sham_fd_order);
#    if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME, (my_crtc () - time1));
#    endif



#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* A operating on psi stored in work2 */
        app_cil_driver (tmp_psiR, work2R, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, sp->hxgrid,
                       sp->hygrid, sp->hzgrid, ct.kohn_sham_fd_order);
        app_cil_driver (tmp_psiI, work2I, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, sp->hxgrid,
                       sp->hygrid, sp->hzgrid, ct.kohn_sham_fd_order);

#    if MD_TIMERS
        rmg_timings (DIAG_APPCIL_TIME, (my_crtc () - time1));
#    endif

        for (idx = 0; idx < get_P0_BASIS(); idx++)
        {
            work1R[idx] = 0.5 * ct.vel * (work1R[idx] - work2R[idx]);
            work1I[idx] = 0.5 * ct.vel * (work1I[idx] - work2I[idx]);
        }




        work1R += get_P0_BASIS();
        work1I += get_P0_BASIS();
	s_psiR += get_P0_BASIS();
	s_psiI += get_P0_BASIS();
    }

    my_free (work2R);

}                               /* subdiag_app_A */






/*Applies B operator to all wavefunctions*/
/*On input b_psi contains s-operator applied to wavefunction*/
void subdiag_app_B (STATE * states, rmg_double_t * b_psiR, rmg_double_t * b_psiI)
{
    int istate, ione=1;
    rmg_double_t *work2R, *work2I, *work1R, *work1I;
#    if MD_TIMERS
    rmg_double_t time1;
#    endif



    my_malloc (work2R, 2 * states[0].sbasis, rmg_double_t);
    work2I = work2R + states[0].sbasis;

    work1R = b_psiR;
    work1I = b_psiI;


    for (istate = 0; istate < ct.num_states; istate++)
    {

        /*Pack S|psi> into smoothing array */
        //pack_ptos (sg_psiR, work1R, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);
        //pack_ptos (sg_psiI, work1I, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID);
	scopy (&states[0].sbasis, work1R, &ione, work2R, &ione);
	scopy (&states[0].sbasis, work1I, &ione, work2I, &ione);


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /*B operating on S|psi> and store in work3 */
        app_cir_driver (work2R, work1R, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, ct.kohn_sham_fd_order);
        app_cir_driver (work2I, work1I, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID, ct.kohn_sham_fd_order);
#    if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME2, (my_crtc () - time1));
#    endif


        work1R += get_P0_BASIS();
        work1I += get_P0_BASIS();
    }

    my_free (work2R);

}                               /* subdiag_app_B */
#endif
#endif
