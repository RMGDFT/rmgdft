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
double *work1;                    /* Smoothing grids */

static void get_qnm_res(double *work_theta, double *kbpsi, double *kbpsi_res);
static void get_nonortho_res(STATE *, double *, STATE *);
extern int it_scf;


void mg_eig(STATE * states, STATE * states1, double *vxc, double *vh,
            double *vnuc, double *rho, double *rhoc, double * vxc_old, double * vh_old)
{
    int istate, ione = 1;
    double t1;
    int st1, ixx, iyy, izz;
    int numst = ct.num_states;
    int item;
    int mxllda2;


    int *ipiv, info;


    int IA =1, JA =1, IB =1, JB=1;

    numst = ct.num_states;

    mxllda2 = MXLLDA * MXLCOL;
    ct.meanres = 0.;
    ct.minres = 100000.0;
    ct.maxres = 0.;

    /* Grab some memory */
    work1 = work_memory;

    void *RT = BeginRmgTimer("3-mg_eig");



    void *RT1a = BeginRmgTimer("3-mg_eig: distribute");
    distribute_to_global(vtot_c, vtot_global);
    EndRmgTimer(RT1a);


   Cpdgemr2d(numst, numst, Hij_00, ione, ione, pct.descb, Hij, ione, ione,
            pct.desca, pct.desca[1]);
    Cpdgemr2d(numst, numst, Bij_00, ione, ione, pct.descb, matB, ione, ione,
            pct.desca, pct.desca[1]);


    void *RT1 = BeginRmgTimer("3-mg_eig: invert MatB");

    my_malloc(ipiv, numst, int);
    /* Compute matrix theta = matB * Hij  */
    pdgetrf(&numst, &numst, matB, &IA, &JA, pct.desca, ipiv, &info);
    if(info !=0)
    { 
        printf("\n error in pdgetrf in mg_eig.c INFO = %d\n", info);
        fflush(NULL);
        exit(0);
    }

    dcopy(&mxllda2, Hij, &ione, theta, &ione);

    pdgetrs("N", &numst, &numst, matB, &IA, &JA, pct.desca, ipiv, 
            theta, &IA, &JA, pct.desca, &info);
    //get_invmat(matB);

    my_free(ipiv);
    EndRmgTimer(RT1);


    void *RT11 = BeginRmgTimer("3-mg_eig: scale theta");
    t1 = 2.0;
    dscal(&mxllda2, &t1, theta, &ione);
    EndRmgTimer(RT11);

    void *RT2 = BeginRmgTimer("3-mg_eig: cpdgemr2d");
    Cpdgemr2d(numst, numst, theta, IA, JA, pct.desca, work_matrix_row, IB, JB,
            pct.descb, pct.desca[1]);
    EndRmgTimer(RT2);




    /* calculate  Theta * S * |states[].psiR > and stored in  states1[].psiR 
     *  sum_j S(phi)_j*(Theta)_ji 
     */

    /*begin shuchun wang */
    void *RT12 = BeginRmgTimer("3-mg_eig: dcopya");
    dcopy(&pct.psi_size, states[ct.state_begin].psiR, &ione,
            states_tem[ct.state_begin].psiR, &ione);
    EndRmgTimer(RT12);

    /*  add q_n,m|beta_n><beta_m|psi> on states_res.psiR */



    void *RT3 = BeginRmgTimer("3-mg_eig: nonortho");
    get_nonortho_res(states, work_matrix_row, states1);
    EndRmgTimer(RT3);
    void *RT4 = BeginRmgTimer("3-mg_eig: qnm");
    get_qnm_res(work_matrix_row, kbpsi, kbpsi_res);

    my_barrier();
    EndRmgTimer(RT4);
    /* end shuchun wang */




    /* Loop over states istate (correction of state istate) */
    void *RT2a = BeginRmgTimer("3-mg_eig: mask");

    for (istate = ct.state_begin; istate < ct.state_end; istate++)
    {
        app_mask(istate, states1[istate].psiR, 0);
    }

    EndRmgTimer(RT2a);

    /* Loop over states istate to compute the whole matrix Hij 
       and all the vectors H|psi> */
    /* calculate the H |phi> on this processor and stored in states1[].psiR[] */

    void *RTa = BeginRmgTimer("3-mg_eig: Hpsi");
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        ixx = states[st1].ixmax - states[st1].ixmin + 1;
        iyy = states[st1].iymax - states[st1].iymin + 1;
        izz = states[st1].izmax - states[st1].izmin + 1;

        /* Generate 2*V*psi and store it  in orbit_tem */
        genvlocpsi(states[st1].psiR, st1, orbit_tem, vtot_global, states);
        t1 = -1.0;
        daxpy(&states[st1].size, &t1, orbit_tem, &ione, states1[st1].psiR, &ione);

        /* Pack psi into smoothing array */
        /*		pack_ptos(sg_orbit, states[st1].psiR, ixx, iyy, izz); 
         */
        /* A operating on psi stored in work2 */
        /*		app_cil(sg_orbit, orbit_tem, ixx, iyy, izz, get_hxgrid(), get_hygrid(), get_hzgrid()); 
         */
        app10_del2(states[st1].psiR, orbit_tem, ixx, iyy, izz, get_hxgrid(), get_hygrid(), get_hzgrid());

        t1 = 1.0;
        daxpy(&states[st1].size, &t1, orbit_tem, &ione, states1[st1].psiR, &ione);


        /*                                                                     
         * Add the contribution of the non-local potential to the 
         * residual 
         */
        /* Get the non-local operator acting on psi and store in nvtot */
        item = (st1 - ct.state_begin) * pct.n_ion_center * ct.max_nl;
    void *RT5 = BeginRmgTimer("3-mg_eig: dnm");
        get_dnmpsi(&states[st1], &kbpsi[item], &kbpsi_res[item], orbit_tem);       /*shuchun wang */
    EndRmgTimer(RT5);


        t1 = 2.0;
        daxpy(&states[st1].size, &t1, orbit_tem, &ione, states1[st1].psiR, &ione);
        app_mask(st1, states1[st1].psiR, 0);
    }                           /* end for st1 = .. */

    EndRmgTimer(RTa);
     

    /*
     *    precond(states1[ct.state_begin].psiR);
     *    t1 = 0.1; 
     *    daxpy(&pct.psi_size, &t1, states1[ct.state_begin].psiR, &ione, states[ct.state_begin].psiR, &ione); 
     */


    /*  SD, Pulay or KAIN method for Linear and Nonlinear equations */

    void *RT6a = BeginRmgTimer("3-mg_eig: scale");
    for (istate = ct.state_begin; istate < ct.state_end; istate++)
    {
        t1 = -1.0;
        dscal(&states1[istate].size, &t1, states1[istate].psiR, &ione);
    }

    EndRmgTimer(RT6a);
    if (ct.restart_mix == 1 || ct.move_centers_at_this_step == 1)
    {
        mix_steps = 0;
        if (pct.gridpe == 0)
            printf("\n restart the orbital mixing at step %d \n", ct.scf_steps);
    }


    void *RT6 = BeginRmgTimer("3-mg_eig: mixing+precond");

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

    EndRmgTimer(RT6);

    void *RT7 = BeginRmgTimer("3-mg_eig: ortho_norm");
    ortho_norm_local(states); 

    EndRmgTimer(RT7);

   // normalize_orbits(states);




    firstflag++;

#if     DEBUG
    print_status(states, vh, vxc, vnuc, vh_old, "before leaving mg_eig.c  ");
    print_state_sum(states1);
#endif
    EndRmgTimer(RT);

}                               /* end mg_eig */

/*-----------------------------------------------------------------*/

static void get_nonortho_res(STATE * states, double *work_theta, STATE * states1)
{
    int i;
    int idx, st1, st2;
    double theta_ion;
    int loop, state_per_proc;
    int num_recv;
    double temp;
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

    int ip1, st1, st2;
    MPI_Status mstatus;
    int ion1, ion2, ion1_global, ion2_global;
    int iip1, iip2;
    int size, proc, proc1, proc2, idx;
    int st11;
    int size_projector;

    size_projector = max_ion_nonlocal * ct.max_nl;
    size = ct.state_per_proc * size_projector;
    for(st1 = 0; st1 < size; st1++) kbpsi_res[st1] = 0.0;

    /* Loop over states on this proce onle 
       (distribution of work AND Aij contributions) */
    proc = pct.gridpe;

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

    void *RT1 = BeginRmgTimer("3-mg_eig: qnm: comm_loop");

    for (idx = 0; idx < kbpsi_num_loop; idx++)
    {


    //    MPI_Barrier(pct.grid_comm);

        proc1 = kbpsi_comm_send[idx];
        proc2 = kbpsi_comm_recv[idx];
        
        MPI_Request request;
        if(proc1 >=0) MPI_Isend(kbpsi, size, MPI_DOUBLE, proc1, idx,
                pct.grid_comm, &request);

        if(proc2 >=0)
            MPI_Recv(kbpsi_comm, size, MPI_DOUBLE, proc2, idx,
                    pct.grid_comm, &mstatus);
        if(proc1 >=0) MPI_Wait(&request, &mstatus);
        if(proc2 < 0) continue;




        void *RT2 = BeginRmgTimer("3-mg_eig: qnm: comm_loop: calcu");

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
        EndRmgTimer(RT2);
    } 
    EndRmgTimer(RT1);

    //    for(st1 = 0; st1 <100; st1++) printf("\n kbpsi_res %d %f", st1, kbpsi_res[st1]);
}


