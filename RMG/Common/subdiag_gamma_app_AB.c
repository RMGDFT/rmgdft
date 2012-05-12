

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

#if !HYBRID_MODEL


#  if GAMMA_PT
/*Applies A operator to all wavefunctions*/
void subdiag_app_A (STATE * states, REAL * a_psi, REAL * s_psi, REAL * vtot_eig)
{
    int kidx, idx, istate, sbasis;
    REAL *sg_twovpsi, *tmp_psi, *work2, *work1;
    STATE *sp;
#    if MD_TIMERS
    REAL time1;
#    endif

    sbasis = states[0].sbasis;
    my_malloc (work2, 2 * sbasis, REAL);
    sg_twovpsi = work2 + sbasis;
    kidx = 0;

    work1 = a_psi;


    for (istate = 0; istate < ct.num_states; istate++)
    {

        sp = &states[istate];
        tmp_psi = sp->psiR;


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* Apply non-local operator to psi and store in work2 */
        app_nls (tmp_psi, NULL, work2, NULL, s_psi, NULL, pct.newsintR_local, NULL, sp->istate, kidx);
#    if MD_TIMERS
        rmg_timings (DIAG_NL_TIME, (my_crtc () - time1));
#    endif



#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* Generate 2*V*psi and store it in a smoothing grid and store in sg_twovpsi */
        genvpsi (tmp_psi, sg_twovpsi, vtot_eig, work2, NULL, 0.0, PX0_GRID, PY0_GRID, PZ0_GRID);
#    if MD_TIMERS
        rmg_timings (DIAG_GENVPSI_TIME, (my_crtc () - time1));
#    endif


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* B operating on 2*V*psi stored in work1 */
        app_cir_driver (sg_twovpsi, work1, PX0_GRID, PY0_GRID, PZ0_GRID, ct.kohn_sham_fd_order);
#    if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME, (my_crtc () - time1));
#    endif


        /* Pack psi into smoothing array */
        //pack_ptos (sg_psi, tmp_psi, PX0_GRID, PY0_GRID, PZ0_GRID);


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* A operating on psi stored in work2 */
        app_cil_driver (tmp_psi, work2, PX0_GRID, PY0_GRID, PZ0_GRID, sp->hxgrid,
                       sp->hygrid, sp->hzgrid, ct.kohn_sham_fd_order);

#    if MD_TIMERS
        rmg_timings (DIAG_APPCIL_TIME, (my_crtc () - time1));
#    endif

        for (idx = 0; idx < P0_BASIS; idx++)
            work1[idx] = 0.5 * ct.vel * (work1[idx] - work2[idx]);




        work1 += P0_BASIS;
        s_psi += P0_BASIS;
    }

    my_free (work2);

}                               /* subdiag_app_A */


/*Applies B operator to all wavefunctions*/
/*On input b_psi contains s-operator applied to wavefunction*/
void subdiag_app_B (STATE * states, REAL * b_psi)
{
    int istate, pbasis, ione=1;
    REAL *work2, *work1;
#    if MD_TIMERS
    REAL time1;
#    endif

    pbasis = states[0].pbasis;


    my_malloc (work2, pbasis, REAL);

    work1 = b_psi;


    for (istate = 0; istate < ct.num_states; istate++)
    {

        /*Pack S|psi> into smoothing array */
        //pack_ptos (sg_psi, work1, PX0_GRID, PY0_GRID, PZ0_GRID);
	scopy (&pbasis, work1, &ione, work2, &ione);


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /*B operating on S|psi> and store in work3 */
        app_cir_driver (work2, work1, PX0_GRID, PY0_GRID, PZ0_GRID, ct.kohn_sham_fd_order);
#    if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME2, (my_crtc () - time1));
#    endif


        work1 += P0_BASIS;
    }

    my_free (work2);

}                               /* subdiag_app_B */


#else

/*Applies A operator to all wavefunctions*/
void subdiag_app_A (STATE * states, REAL * a_psiR, REAL * a_psiI, REAL * s_psiR, REAL * s_psiI, REAL * vtot_eig)
{
    int kidx, idx, istate, sbasis;
    REAL *sg_twovpsiR, *sg_twovpsiI, *tmp_psiR, *tmp_psiI, *work2R, *work2I,
        *work1R, *work1I;
    REAL *gx, *gy, *gz, *kdr;
    STATE *sp;
#    if MD_TIMERS
    REAL time1;
#    endif



    sbasis = states[0].sbasis;
    my_malloc (work2R, 8 * sbasis, REAL);
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
        /*pack_ptos (sg_psiR, tmp_psiR, PX0_GRID, PY0_GRID, PZ0_GRID);
        pack_ptos (sg_psiI, tmp_psiI, PX0_GRID, PY0_GRID, PZ0_GRID);*/


        /* Apply the gradient operator to psi */
        app_grad (tmp_psiI, (P0_GRID *) gx, (P0_GRID *) gy, (P0_GRID *) gz);
        for (idx = 0; idx < P0_BASIS; idx++)
            kdr[idx] = (ct.kp[sp->kidx].kvec[0] * gx[idx] +
                        ct.kp[sp->kidx].kvec[1] * gy[idx] + ct.kp[sp->kidx].kvec[2] * gz[idx]);




#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* Generate 2*V*psi and store in sg_twovpsi */
        genvpsi (tmp_psiR, sg_twovpsiR, vtot_eig, work2R, kdr, ct.kp[sp->kidx].kmag, PX0_GRID,
                 PY0_GRID, PZ0_GRID);
#    if MD_TIMERS
        rmg_timings (DIAG_GENVPSI_TIME, (my_crtc () - time1));
#    endif


        /* Apply the gradient operator to psi */
        app_grad (tmp_psiR, (P0_GRID *) gx, (P0_GRID *) gy, (P0_GRID *) gz);
        for (idx = 0; idx < P0_BASIS; idx++)
            kdr[idx] = -(ct.kp[sp->kidx].kvec[0] * gx[idx] +
                         ct.kp[sp->kidx].kvec[1] * gy[idx] + ct.kp[sp->kidx].kvec[2] * gz[idx]);


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* Generate 2 * V * psiI and store it in a smoothing grid and store in sg_twovpsiI */
        genvpsi (tmp_psiI, sg_twovpsiI, vtot_eig, work2I, kdr, ct.kp[sp->kidx].kmag, PX0_GRID,
                 PY0_GRID, PZ0_GRID);
#    if MD_TIMERS
        rmg_timings (DIAG_GENVPSI_TIME, (my_crtc () - time1));
#    endif




#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* B operating on 2*V*psi stored in work1 */
        app_cir_driver (sg_twovpsiR, work1R, PX0_GRID, PY0_GRID, PZ0_GRID, ct.kohn_sham_fd_order);
        app_cir_driver (sg_twovpsiI, work1I, PX0_GRID, PY0_GRID, PZ0_GRID, ct.kohn_sham_fd_order);
#    if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME, (my_crtc () - time1));
#    endif



#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* A operating on psi stored in work2 */
        app_cil_driver (tmp_psiR, work2R, PX0_GRID, PY0_GRID, PZ0_GRID, sp->hxgrid,
                       sp->hygrid, sp->hzgrid, ct.kohn_sham_fd_order);
        app_cil_driver (tmp_psiI, work2I, PX0_GRID, PY0_GRID, PZ0_GRID, sp->hxgrid,
                       sp->hygrid, sp->hzgrid, ct.kohn_sham_fd_order);

#    if MD_TIMERS
        rmg_timings (DIAG_APPCIL_TIME, (my_crtc () - time1));
#    endif

        for (idx = 0; idx < P0_BASIS; idx++)
        {
            work1R[idx] = 0.5 * ct.vel * (work1R[idx] - work2R[idx]);
            work1I[idx] = 0.5 * ct.vel * (work1I[idx] - work2I[idx]);
        }




        work1R += P0_BASIS;
        work1I += P0_BASIS;
	s_psiR += P0_BASIS;
	s_psiI += P0_BASIS;
    }

    my_free (work2R);

}                               /* subdiag_app_A */






/*Applies B operator to all wavefunctions*/
/*On input b_psi contains s-operator applied to wavefunction*/
void subdiag_app_B (STATE * states, REAL * b_psiR, REAL * b_psiI)
{
    int istate, ione=1;
    REAL *work2R, *work2I, *work1R, *work1I;
#    if MD_TIMERS
    REAL time1;
#    endif



    my_malloc (work2R, 2 * states[0].sbasis, REAL);
    work2I = work2R + states[0].sbasis;

    work1R = b_psiR;
    work1I = b_psiI;


    for (istate = 0; istate < ct.num_states; istate++)
    {

        /*Pack S|psi> into smoothing array */
        //pack_ptos (sg_psiR, work1R, PX0_GRID, PY0_GRID, PZ0_GRID);
        //pack_ptos (sg_psiI, work1I, PX0_GRID, PY0_GRID, PZ0_GRID);
	scopy (&states[0].sbasis, work1R, &ione, work2R, &ione);
	scopy (&states[0].sbasis, work1I, &ione, work2I, &ione);


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /*B operating on S|psi> and store in work3 */
        app_cir_driver (work2R, work1R, PX0_GRID, PY0_GRID, PZ0_GRID, ct.kohn_sham_fd_order);
        app_cir_driver (work2I, work1I, PX0_GRID, PY0_GRID, PZ0_GRID, ct.kohn_sham_fd_order);
#    if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME2, (my_crtc () - time1));
#    endif


        work1R += P0_BASIS;
        work1I += P0_BASIS;
    }

    my_free (work2R);

}                               /* subdiag_app_B */
#endif
#endif
