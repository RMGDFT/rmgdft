/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
                            mg_eig.c 
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "my_scalapack.h"

/* Flag for projector mixing */
int firstflag = FALSE;
static int mix_steps;
rmg_double_t *work1;                    /* Smoothing grids */

static void get_qnm_res(double *work_theta, double *kbpsi, double *kbpsi_res);
static void get_nonortho_res(STATE *, double *, STATE *);
extern int it_scf;


void mg_eig(STATE * states, STATE * states1, double *vxc, double *vh,
            double *vnuc, double *rho, double *rhoc, rmg_double_t * vxc_old, rmg_double_t * vh_old)
{
    int idx, istate, ione = 1;
    double diag, t1, d1, gamma;
    STATE *sp, *sp1;
    int st1, ixx, iyy, izz;
    char side = 'l', uplo = 'l';
    int n2 = ct.num_states * ct.num_states, numst = ct.num_states;
    int item;
    rmg_double_t alpha1, alpha2;
    rmg_double_t tem;
    int mxllda2;
    rmg_double_t A11, A12, A21, A22, c1, c2, b1, b2;
    rmg_double_t x1r1, x2r1, x3r1, x1r2, x2r2, x3r2, x1r3, x2r3, x3r3;


    rmg_double_t tem_luw = -1.0;


    int IA =1, JA =1, IB =1, JB=1;

    numst = ct.num_states;

    mxllda2 = MXLLDA * MXLCOL;
    diag = -1. / ct.Ac;
    ct.meanres = 0.;
    ct.minres = 100000.0;
    ct.maxres = 0.;

    /* Grab some memory */
    work1 = work_memory;



    distribute_to_global(vtot_c, vtot_global);


    get_invmat(matB);

    /* Compute matrix theta = matB * Hij  */
    dsymm_dis(&side, &uplo, &numst, matB, Hij, theta);

    t1 = 2.0;
    sscal(&mxllda2, &t1, theta, &ione);

    Cpdgemr2d(numst, numst, theta, IA, JA, pct.desca, work_matrix_row, IB, JB,
            pct.descb, pct.desca[1]);




    /* calculate  Theta * S * |states[].psiR > and stored in  states1[].psiR 
     *  sum_j S(phi)_j*(Theta)_ji 
     */

    /*begin shuchun wang */
    scopy(&pct.psi_size, states[ct.state_begin].psiR, &ione,
            states_tem[ct.state_begin].psiR, &ione);

    /*  add q_n,m|beta_n><beta_m|psi> on states_res.psiR */



    get_nonortho_res(states, work_matrix_row, states1);
    get_qnm_res(work_matrix_row, kbpsi, kbpsi_res);

    my_barrier();
    /* end shuchun wang */




    /* Loop over states istate (correction of state istate) */

    for (istate = ct.state_begin; istate < ct.state_end; istate++)
    {
        app_mask(istate, states1[istate].psiR, 0);
    }


    /* Loop over states istate to compute the whole matrix Hij 
       and all the vectors H|psi> */
    /* calculate the H |phi> on this processor and stored in states1[].psiR[] */

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        ixx = states[st1].ixmax - states[st1].ixmin + 1;
        iyy = states[st1].iymax - states[st1].iymin + 1;
        izz = states[st1].izmax - states[st1].izmin + 1;

        /* Generate 2*V*psi and store it  in orbit_tem */
        genvlocpsi(states[st1].psiR, st1, orbit_tem, vtot_global, states);
        t1 = -1.0;
        saxpy(&states[st1].size, &t1, orbit_tem, &ione, states1[st1].psiR, &ione);

        /* Pack psi into smoothing array */
        /*		pack_ptos(sg_orbit, states[st1].psiR, ixx, iyy, izz); 
         */
        /* A operating on psi stored in work2 */
        /*		app_cil(sg_orbit, orbit_tem, ixx, iyy, izz, get_hxgrid(), get_hygrid(), get_hzgrid()); 
         */
        app10_del2(states[st1].psiR, orbit_tem, ixx, iyy, izz, get_hxgrid(), get_hygrid(), get_hzgrid());

        t1 = 1.0;
        saxpy(&states[st1].size, &t1, orbit_tem, &ione, states1[st1].psiR, &ione);


        /*                                                                     
         * Add the contribution of the non-local potential to the 
         * residual 
         */
        /* Get the non-local operator acting on psi and store in nvtot */
        item = (st1 - ct.state_begin) * pct.n_ion_center * ct.max_nl;
        get_dnmpsi(&states[st1], &kbpsi[item], &kbpsi_res[item], orbit_tem);       /*shuchun wang */


        t1 = 2.0;
        saxpy(&states[st1].size, &t1, orbit_tem, &ione, states1[st1].psiR, &ione);
        app_mask(st1, states1[st1].psiR, 0);
    }                           /* end for st1 = .. */

     

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
        if (pct.gridpe == 0)
            printf("\n restart the orbital mixing at step %d \n", ct.scf_steps);
    }



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


    mix_steps++;

     	ortho_norm_local(states); 

   // normalize_orbits(states);




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
    int idx, st1, st2;
    rmg_double_t theta_ion;
    rmg_double_t *psi_pointer,*psi3,*psi2, *psi1;
    int loop, state_per_proc, proc1, proc2;
    int num_send, num_recv, num_sendrecv, size1, size2;
    rmg_double_t temp;
    MPI_Status mstatus;
    MPI_Request  mr_send,*mr_recv;
    int st11;


    state_per_proc = ct.state_per_proc + 2;

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        for (idx = 0; idx < states[st1].size; idx++)
            states1[st1].psiR[idx] = 0.0;

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        st11 = st1-ct.state_begin;
        for (st2 = ct.state_begin; st2 < ct.state_end; st2++)
            if (state_overlap_or_not[st11 * ct.num_states + st2] == 1)
            {
                temp = work_theta[st11 * ct.num_states + st2];
                theta_phi_new(st1, st2, temp, states[st2].psiR, states1[st1].psiR, 0, states);
            }
    }

    /*  send overlaped orbitals to other PEs */
    my_barrier();




    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {

        num_recv = recv_from1[loop * state_per_proc + 1];

        for (i = 0; i < num_recv; i++)
        {
            st2 = recv_from1[loop * state_per_proc + i + 2];
            for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
            {
                st11 = st1-ct.state_begin;
                if (state_overlap_or_not[st11 * ct.num_states + st2] == 1)
                {
                    theta_ion = work_theta[st11 * ct.num_states + st2];
                    theta_phi_new(st1, st2, theta_ion, states[st2].psiR, states1[st1].psiR, 0, states);
                }
            }
        }

    }

}

void get_qnm_res(double *work_theta, double *kbpsi, double *kbpsi_res)
{

    int ion, ip1, ip2, st1, st2, ist;
    MPI_Status mstatus;
    int ion1, ion2, ion1_global, ion2_global;
    int iip1, iip2, iip1a, iip2a;
    int size, proc, proc1, proc2, idx;
    int nh;
    int st11;
    double one = 1.0, zero= 0.0;
    int size_projector, numst_thispe;

    size_projector = max_ion_nonlocal * ct.max_nl;
    size = ct.state_per_proc * size_projector;
    for(st1 = 0; st1 < size; st1++) kbpsi_res[st1] = 0.0;

    /* Loop over states on this proce onle 
       (distribution of work AND Aij contributions) */
    proc = pct.gridpe;
    numst_thispe = ct.state_end - ct.state_begin;

    for (ion1 = 0; ion1 < num_nonlocal_ion[proc]; ion1++)
        for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        {
            st11 = st1 - ct.state_begin;

            iip1 = (st1 - state_begin[proc]) * num_nonlocal_ion[proc] * ct.max_nl;
            iip1 += ion1 * ct.max_nl;

            for (st2 = state_begin[proc]; st2 < state_end[proc]; st2++)
            {
                //if(state_overlap_or_not[st11 * ct.num_states + st2] == 1)
                {
                    iip2 = (st2 - state_begin[proc]) * num_nonlocal_ion[proc] * ct.max_nl;
                    iip2 += ion1 * ct.max_nl;
                    for (ip1 = 0; ip1 < ct.max_nl; ip1++)
                        kbpsi_res[iip1 + ip1] += kbpsi[iip2 + ip1] * work_theta[st11 * ct.num_states + st2];
                }
            }

        }         


    for (idx = 1; idx < NPES; idx++)
    {

        proc1 = pct.gridpe + idx;
        if (proc1 >= NPES)
            proc1 = proc1 - NPES;
        proc2 = pct.gridpe - idx;
        if (proc2 < 0)
            proc2 += NPES;


        MPI_Sendrecv(kbpsi, size, MPI_DOUBLE, proc1, idx, kbpsi_comm, size,
                MPI_DOUBLE, proc2, idx, pct.grid_comm, &mstatus);



        for (ion1 = 0; ion1 < num_nonlocal_ion[proc]; ion1++)
            for (ion2 = 0; ion2 < num_nonlocal_ion[proc2]; ion2++)
            {
                ion1_global = ionidx_allproc[proc * max_ion_nonlocal + ion1];
                ion2_global = ionidx_allproc[proc2 * max_ion_nonlocal + ion2];

                if (ion1_global == ion2_global)
                {
                    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                    {
                        st11 = st1 - ct.state_begin;

                        iip1 = (st1 - state_begin[proc]) * num_nonlocal_ion[proc] * ct.max_nl;
                        iip1 += ion1 * ct.max_nl;

                        for (st2 = state_begin[proc2]; st2 < state_end[proc2]; st2++)
                        {

                            iip2 = (st2 - state_begin[proc2]) * num_nonlocal_ion[proc2] * ct.max_nl;
                            iip2 += ion2 * ct.max_nl;
                            for (ip1 = 0; ip1 < ct.max_nl; ip1++)
                                kbpsi_res[iip1+ip1] += kbpsi_comm[iip2+ip1] * work_theta[st11 * ct.num_states + st2];
                        }

                    }         
                }            

            }
    } 
}
