/************************** SVN Revision Information **************************
 **    $Id: mg_eig.c 1108 2010-03-30 21:23:14Z luw $    **
******************************************************************************/
 
/*
                            mg_eig.c 
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"

#include "my_scalapack.h"

/* Flag for projector mixing */
int firstflag = FALSE;
static int mix_steps;
REAL *work1;                    /* Smoothing grids */

static void get_nonortho_res(STATE *, double *, STATE *);
extern int it_scf;


void mg_eig(STATE * states, STATE * states1, double *vxc, double *vh,
            double *vnuc, double *rho, double *rhoc, REAL * vxc_old, REAL * vh_old)
{
    int idx, istate, ione = 1;
    double diag, t1, d1, gamma;
    STATE *sp, *sp1;
    double time1, time2;
    int st1, ixx, iyy, izz;
    char side = 'l', uplo = 'l';
    int n2 = ct.num_states * ct.num_states, numst = ct.num_states;
    int item;
    REAL alpha1, alpha2;
    REAL tem;
    int mxllda2;
    REAL A11, A12, A21, A22, c1, c2, b1, b2;
    REAL x1r1, x2r1, x3r1, x1r2, x2r2, x3r2, x1r3, x2r3, x3r3;


    REAL tem_luw = -1.0;


    int IA =1, JA =1, IB =1, JB=1;

    numst = ct.num_states;

    mxllda2 = MXLLDA * MXLLDA;
    diag = -1. / ct.Ac;
    ct.meanres = 0.;
    ct.minres = 100000.0;
    ct.maxres = 0.;

    /* Grab some memory */
    work1 = work_memory;


    get_ddd(vtot);

    time1 = my_crtc();
    get_vtot_psi(vtot_c, vtot);

    distribute_to_global(vtot_c, vtot_global);
    time2 = my_crtc();
    rmg_timings(POTFC_TIME, time2 - time1, 0);


    gamma = get_gamma(vtot_c, states[0].eig);

    if (pct.thispe == 0)
    {
        printf("\n time step for low frequencies corrections = %e\n", gamma);
        printf(" levels= %d\n", ct.eig_parm.levels);
    }

    get_invmat(matB);

    /* Compute matrix theta = matB * Hij  */
    time1 = my_crtc();
    dsymm_dis(&side, &uplo, &numst, matB, Hij, theta);

    t1 = 2.0;
    sscal(&mxllda2, &t1, theta, &ione);

    Cpdgemr2d(numst, numst, theta, IA, JA, pct.desca, work_matrix_row, IB, JB,
            pct.descb, pct.descb[1]);


    time2 = my_crtc();
    rmg_timings(THETA_TIME, time2 - time1, 0);


    /* calculate  Theta * S * |states[].psiR > and stored in  states1[].psiR 
     *  sum_j S(phi)_j*(Theta)_ji 
     */

    /*begin shuchun wang */
    scopy(&pct.psi_size, states[ct.state_begin].psiR, &ione,
            states_tem[ct.state_begin].psiR, &ione);

    /*  add q_n,m|beta_n><beta_m|psi> on states_res.psiR */


    time1 = my_crtc();

    for (istate = ct.state_begin; istate < ct.state_end; istate++)
    {

        /* calculate the q_n,m |beta_n><beta_m|psi>  on this processor and stored in states1[].psiR[] */
        item = (istate - ct.state_begin) * pct.n_ion_center * ct.max_nl;
        get_qnmpsi(&states[istate], &kbpsi[item], orbit_tem);
        t1 = 1.0;
        saxpy(&states[istate].size, &t1, orbit_tem, &ione, states_tem[istate].psiR, &ione);
    }

    time2 = my_crtc();
    rmg_timings(QNMPSI_TIME, time2 - time1, 0);


    get_nonortho_res(states_tem, work_matrix_row, states1);
    my_barrier();
    time1 = my_crtc();
    rmg_timings(NONRES_TIME, time1 - time2, 0);
    /* end shuchun wang */




    /* Loop over states istate (correction of state istate) */

    for (istate = ct.state_begin; istate < ct.state_end; istate++)
    {
        ixx = states[istate].ixmax - states[istate].ixmin + 1;
        iyy = states[istate].iymax - states[istate].iymin + 1;
        izz = states[istate].izmax - states[istate].izmin + 1;

        app_mask(istate, states1[istate].psiR, 0);
        pack_ptos(sg_orbit, states1[istate].psiR, ixx, iyy, izz);
        fill_orbit_borders(sg_orbit, ixx, iyy, izz);

        /* res = B*sg_res */
        app_cir_0(sg_orbit, orbit_tem, ixx, iyy, izz);
        app_mask(istate, orbit_tem, 0);

        for (idx = 0; idx < ixx * iyy * izz; idx++)
            states1[istate].psiR[idx] = orbit_tem[idx];

        app_mask(istate, states1[istate].psiR, 0);
    }


    /* Loop over states istate to compute the whole matrix Hij 
       and all the vectors H|psi> */
    /* calculate the H |phi> on this processor and stored in states1[].psiR[] */

    time1 = my_crtc();
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        ixx = states[st1].ixmax - states[st1].ixmin + 1;
        iyy = states[st1].iymax - states[st1].iymin + 1;
        izz = states[st1].izmax - states[st1].izmin + 1;

        /* Generate 2*V*psi and store it  in orbit_tem */
        genvlocpsi(states[st1].psiR, st1, orbit_tem, vtot_global, states);
        pack_ptos(sg_orbit, orbit_tem, ixx, iyy, izz);
        fill_orbit_borders(sg_orbit, ixx, iyy, izz);

        /* B operating on 2*V*psi stored in work1 */
        app_cir_0(sg_orbit, orbit_tem, ixx, iyy, izz);
        t1 = -1.0;
        saxpy(&states[st1].size, &t1, orbit_tem, &ione, states1[st1].psiR, &ione);

        /* Pack psi into smoothing array */
        /*		pack_ptos(sg_orbit, states[st1].psiR, ixx, iyy, izz); 
         *		fill_orbit_borders(sg_orbit, ixx, iyy, izz); 
         */
        /* A operating on psi stored in work2 */
        /*		app_cil(sg_orbit, orbit_tem, ixx, iyy, izz, ct.hxgrid, ct.hygrid, ct.hzgrid); 
         */
        app10_del2(states[st1].psiR, orbit_tem, ixx, iyy, izz, ct.hxgrid, ct.hygrid, ct.hzgrid);

        for (idx = 0; idx < ixx * iyy * izz; idx++)
            states1[st1].psiR[idx] += orbit_tem[idx];

        app_mask(st1, states1[st1].psiR, 0);

        /*                                                                     
         * Add the contribution of the non-local potential to the 
         * residual 
         */
        /* Get the non-local operator acting on psi and store in nvtot */
        item = (st1 - ct.state_begin) * pct.n_ion_center * ct.max_nl;
        get_dnmpsi(&states[st1], &kbpsi[item], orbit_tem);       /*shuchun wang */

        pack_ptos(sg_orbit, orbit_tem, ixx, iyy, izz);
        fill_orbit_borders(sg_orbit, ixx, iyy, izz);

        /* B operating on 2*VNL*psi stored in orbit_tem */
        app_cir_0(sg_orbit, orbit_tem, ixx, iyy, izz);
        app_mask(st1, orbit_tem, 0);

        t1 = -2.0;
        saxpy(&states[st1].size, &t1, orbit_tem, &ione, states1[st1].psiR, &ione);
    }                           /* end for st1 = .. */
    time2 = my_crtc();
    rmg_timings(HPSI_TIME, time2 - time1, 0);

    /*  print_sum(pct.psi_size, states1[ct.state_begin].psiR, "states1 sum "); 
     */

    /*
     *    precond(states1[ct.state_begin].psiR);
     *    t1 = 0.1; 
     *    saxpy(&pct.psi_size, &t1, states1[ct.state_begin].psiR, &ione, states[ct.state_begin].psiR, &ione); 
     */


    /*  SD, Pulay or KAIN method for Linear and Nonlinear equations */

    for (istate = ct.state_begin; istate < ct.state_end; istate++)
    {
        t1 = -1.0;
        sscal(&states1[istate].size, &t1, states1[istate].psiR, &ione);
    }

    if (ct.restart_mix == 1 || ct.move_centers_at_this_step == 1)
    {
        mix_steps = 0;
        if (pct.thispe == 0)
            printf("\n restart the orbital mixing at step %d \n", ct.scf_steps);
    }


    time1 = my_crtc();

    switch (ct.mg_method)
    {
        case 0:
            sd(ct.scf_steps, pct.psi_size, states[ct.state_begin].psiR, states1[ct.state_begin].psiR);
            break;
        case 1:
            pulay(mix_steps, pct.psi_size, states[ct.state_begin].psiR,
                    states1[ct.state_begin].psiR, ct.mg_steps, 1);
            break;
        case 2:
            kain(mix_steps, pct.psi_size, states[ct.state_begin].psiR,
                    states1[ct.state_begin].psiR, ct.mg_steps);
            break;
        case 3:
            pulay_weighted(mix_steps, pct.psi_size, states[ct.state_begin].psiR,
                    states1[ct.state_begin].psiR, ct.mg_steps, 100, 0.5, 1);
            break;
        default:
            error_handler("Undefined mg_method ");
    }

    time2 = my_crtc();
    rmg_timings(MIXPSI_TIME, time2 - time1, 0);

    mix_steps++;

    /* 	ortho_norm_local(states); 
     */

    normalize_orbits(states);



    time2 = my_crtc();
    d1 = time2 - time1;
    rmg_timings(MG_TIME, d1, 0);

    firstflag++;

#if     DEBUG
    print_status(states, vh, vxc, vnuc, vh_old, "before leaving mg_eig.c  ");
    print_state_sum(states1);
#endif

}                               /* end mg_eig */

/*-----------------------------------------------------------------*/

static void get_nonortho_res(STATE * states, double *work_theta, STATE * states1)
{
    int i,ii,max_ii;
    double time1, time2;
    int idx, st1, st2;
    REAL theta_ion;
    REAL *psi_pointer,*psi3,*psi2, *psi1;
    int loop, state_per_proc, proc1, proc2;
    int num_send, num_recv, num_sendrecv, size1, size2;
    REAL temp;
    MPI_Status mstatus;
    MPI_Request  mr_send,*mr_recv;
    int st11;


    state_per_proc = ct.state_per_proc + 2;
    psi2 = orbit_tem;
    my_malloc_init(psi3, ct.max_orbit_size, REAL);

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        for (idx = 0; idx < states[st1].size; idx++)
            states1[st1].psiR[idx] = 0.0;

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        for (st2 = ct.state_begin; st2 < ct.state_end; st2++)
            if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
            {
                st11 = st1-ct.state_begin;
                temp = work_theta[st11 * ct.num_states + st2];
                theta_phi_new(st1, st2, temp, states[st2].psiR, states1[st1].psiR, 0, states);
            }

    /*  send overlaped orbitals to other PEs */
    my_barrier();



    max_ii = 0;
    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {
        num_send = send_to1[loop * state_per_proc + 1];
        num_recv = recv_from1[loop * state_per_proc + 1];
        max_ii = max(max_ii, num_send);
        max_ii = max(max_ii, num_recv);
    }

    max_ii = int_max_all(max_ii);

    ii = num_sendrecv_loop1 * (max_ii +10) +1;

    my_calloc(mr_recv, ii, MPI_Request);

    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {

        proc1 = send_to1[loop * state_per_proc];
        proc2 = recv_from1[loop * state_per_proc];
        num_send = send_to1[loop * state_per_proc + 1];
        num_recv = recv_from1[loop * state_per_proc + 1];

        ii = loop * max_ii+1;
        if(num_recv>0)
        {
            i=0;
            ii++;
            st2 = recv_from1[loop * state_per_proc + i + 2];
            size2 = states[st2].size;
            if(ii %2 ==0) MPI_Irecv(psi2, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
            if(ii %2 ==1) MPI_Irecv(psi3, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
        }



        ii = loop * max_ii+1;
        for (i = 0; i < num_send; i++)
        {
            ii++;
            st1 = send_to1[loop * state_per_proc + i + 2];
            psi1 = states[st1].psiR;
            size1 = states[st1].size;
            MPI_Isend(psi1, size1, MPI_DOUBLE, proc1, ii, MPI_COMM_WORLD, &mr_send);
            MPI_Request_free(&mr_send);
        }


        ii = loop * max_ii+2;
        for (i = 1; i < num_recv+1; i++)
        {
            ii++;
            MPI_Wait(&mr_recv[ii-1],&mstatus); 
            st2 = recv_from1[loop * state_per_proc + i-1 + 2];
            size2 = states[st2].size;
            if(i != num_recv && ii%2 ==0) MPI_Irecv(psi2, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
            if(i != num_recv && ii%2 ==1) MPI_Irecv(psi3, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
            if(ii%2 ==0) psi_pointer = psi3;
            if(ii%2 ==1) psi_pointer = psi2;
            for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
                {
                    st11 = st1-ct.state_begin;
                    theta_ion = work_theta[st11 * ct.num_states + st2];
                    theta_phi_new(st1, st2, theta_ion, psi_pointer, states1[st1].psiR, 0, states);
                }
        }

    }

    my_free(mr_recv);
    my_free(psi3);
}
