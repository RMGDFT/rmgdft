/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>



#include "make_conf.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"

#include "prototypes_on.h"
#include "init_var.h"
#include "my_scalapack.h"
#include "transition.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"

/* Flag for projector mixing */
static int mix_steps;

static void get_qnm_res(double *work_theta);
static void get_nonortho_res(STATE *, double *, STATE *);
static void get_dnmpsi(STATE *);

void OrbitalOptimize(STATE * states, STATE * states1, double *vxc, double *vh,
            double *vnuc, double *rho, double *rhoc, double * vxc_old, double * vh_old)
{
    int ione = 1;
    double t1;
    int order = ct.kohn_sham_fd_order;
    double hxgrid = Rmg_G->get_hxgrid(1);
    double hygrid = Rmg_G->get_hygrid(1);
    double hzgrid = Rmg_G->get_hzgrid(1);

    int item = (ct.max_orbit_nx+order) *(ct.max_orbit_ny+order) *(ct.max_orbit_nz+order);
    double *orbital_border = new double[2*item];
    double *orbit_tem = new double[2*item];

    FiniteDiff FD(&Rmg_L);

    RmgTimer *RT = new RmgTimer("3-OrbitalOptimize");


    RmgTimer *RT1a = new RmgTimer("3-OrbitalOptimize: distribute");
    distribute_to_global(vtot_c, vtot_global);
    delete(RT1a);


    /* calculate  Theta * S * |states[].psiR > and stored in  states1[].psiR 
     *  sum_j S(phi)_j*(Theta)_ji 
     */

    /*begin shuchun wang */
    RmgTimer *RT12 = new RmgTimer("3-OrbitalOptimize: dcopya");
    dcopy(&pct.psi_size, states[ct.state_begin].psiR, &ione, states_tem[ct.state_begin].psiR, &ione);
    delete(RT12);

    /*  add q_n,m|beta_n><beta_m|psi> on states_res.psiR */
    RmgTimer *RT3 = new RmgTimer("3-OrbitalOptimize: nonortho");
    get_nonortho_res(states, theta, states1);
    delete(RT3);

    RmgTimer *RT4 = new RmgTimer("3-OrbitalOptimize: qnm");
    get_qnm_res(theta);
    delete(RT4);
    /* end shuchun wang */


    /* Loop over states istate to compute the whole matrix Hij 
       and all the vectors H|psi> */
    /* calculate the H |phi> on this processor and stored in states1[].psiR[] */

    RmgTimer *RTa = new RmgTimer("3-OrbitalOptimize: Hpsi");
    for (int st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        STATE *sp = &states[st1];
        STATE *sp1 = &states1[st1];
        int ixx = states[st1].orbit_nx;
        int iyy = states[st1].orbit_ny;
        int izz = states[st1].orbit_nz;

        /* Generate 2*V*psi and store it  in orbit_tem */
        genvlocpsi(sp->psiR, st1, orbit_tem, vtot_global, states);

        t1 = -1.0;
        daxpy(&sp->size, &t1, orbit_tem, &ione, sp1->psiR, &ione);

        /* A operating on psi stored in orbit_tem */
        if(sp->radius > 0)
        {
            FillOrbitalBorders(orbital_border, sp->psiR, ixx, iyy, izz, order);
        }
        else
        {
            Rmg_T->trade_imagesx_central_local(sp->psiR, orbital_border, ixx, iyy, izz, order/2);
        }
        FD.app8_del2 (orbital_border, orbit_tem, ixx, iyy, izz, hxgrid, hygrid, hzgrid);

        t1 = 1.0;
        daxpy(&sp->size, &t1, orbit_tem, &ione, sp1->psiR, &ione);

    }                           /* end for st1 = .. */
    RmgTimer *RT5 = new RmgTimer("3-OrbitalOptimize: dnm");
    get_dnmpsi(states1);
    delete(RT5);


    delete(RTa);


    /*
     *    precond(states1[ct.state_begin].psiR);
     *    t1 = 0.1; 
     *    daxpy(&pct.psi_size, &t1, states1[ct.state_begin].psiR, &ione, states[ct.state_begin].psiR, &ione); 
     */


    /*  SD, Pulay or KAIN method for Linear and Nonlinear equations */

    RmgTimer *RT6a = new RmgTimer("3-OrbitalOptimize: scale");
    for (int istate = ct.state_begin; istate < ct.state_end; istate++)
    {
        t1 = -1.0;
        dscal(&states1[istate].size, &t1, states1[istate].psiR, &ione);
    }

    delete(RT6a);
    if (ct.restart_mix == 1 || ct.move_centers_at_this_step == 1)
    {
        mix_steps = 0;
        if (pct.gridpe == 0)
            printf("\n restart the orbital mixing at step %d \n", ct.scf_steps);
    }


    RmgTimer *RT6 = new RmgTimer("3-OrbitalOptimize: mixing+precond");

    double gamma = -0.5;
    switch (ct.mg_method)
    {
        case 0:
            Precond(states1[ct.state_begin].psiR);
            daxpy(&pct.psi_size, &gamma, states1[ct.state_begin].psiR, &ione, states[ct.state_begin].psiR, &ione);
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
            printf("\n Undefined mg_method \n ");
            fflush(NULL);
            exit(0);
    }


    mix_steps++;

    delete(RT6);

    RmgTimer *RT7 = new RmgTimer("3-OrbitalOptimize: ortho_norm");
    //ortho_norm_local(states); 

    delete(RT7);

    normalize_orbits(states);





#if     DEBUG
    print_status(states, vh, vxc, vnuc, vh_old, "before leaving OrbitalOptimize.c  ");
    print_state_sum(states1);
#endif
    delete [] orbit_tem;
    delete [] orbital_border;
    delete(RT);

} 

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

void get_qnm_res(double *work_theta)
{

    int st1, st2;
    int ion1;
    int st11;
    int num_prj, num_orb, tot_orb, idx1, idx2;
    int max_orb;
    double one = 1.0, zero = 0.0, *work_mat;

    max_orb = 0;

    for (ion1 = 0; ion1 < pct.n_ion_center; ion1++)
    {
        tot_orb = Kbpsi_str.orbital_index[ion1].size();
        max_orb = std::max(max_orb, tot_orb);
    }

    work_mat = new double[(ct.state_end-ct.state_begin) *max_orb];


    for (ion1 = 0; ion1 < pct.n_ion_center; ion1++)
    {

        num_prj = pct.prj_per_ion[pct.ionidx[ion1]];
        if (num_prj == 0) continue;
        num_orb = Kbpsi_str.num_orbital_thision[ion1]; 
        tot_orb = Kbpsi_str.orbital_index[ion1].size();

        for(idx1 = 0; idx1 < num_orb; idx1++)
        {
            st1 = Kbpsi_str.orbital_index[ion1][idx1];
            st11 = st1 - ct.state_begin;

            for(idx2 = 0; idx2 < tot_orb; idx2++)
            {
                st2 = Kbpsi_str.orbital_index[ion1][idx2];
                work_mat[idx1 * tot_orb + idx2] = work_theta[st11 *ct.num_states + st2];
            }
        }

        //  set the length of vector and set their value to 0.0
        Kbpsi_str.kbpsi_res_ion[ion1].resize(num_orb * num_prj);

        dgemm("N", "N", &num_prj, &num_orb, &tot_orb,  &one, Kbpsi_str.kbpsi_ion[ion1].data(), &num_prj,
                work_mat, &tot_orb, &zero, Kbpsi_str.kbpsi_res_ion[ion1].data(), &num_prj);


    }         

    delete [] work_mat;


}


void get_dnmpsi(STATE *states1)
{
    int ion;
    double *prjptr;
    int ion2, st0, st1;
    double *ddd, *qnm_weight;
    double *qqq;
    double *prj_sum;

    double one=1.0, zero=0.0, mtwo=-2.0, *work_kbpsi; 
    int ione=1, num_orb, num_prj;

    /*
     *  dnm_weght[i] = sum_j dnm[i, j] *<beta_j|phi> 
     * prj_sum(r) = sum_i nm_weight[i] *beta_i(r)  
     */

    qnm_weight = new double[ct.max_nl];
    prj_sum = new double[ct.max_nlpoints];
    work_kbpsi = new double[ct.max_nl * (ct.state_end-ct.state_begin)];


    prjptr = projectors;

    for (ion2 = 0; ion2 < pct.n_ion_center; ion2++)
    {
        ion = pct.ionidx[ion2];
        num_prj = pct.prj_per_ion[ion];
        if (num_prj == 0) continue;
        qqq = pct.qqq[ion];
        ddd = pct.dnmI[ion];
        num_orb = Kbpsi_str.num_orbital_thision[ion2]; 

        dgemm("N", "N", &num_prj, &num_orb, &num_prj, &one, qqq, &num_prj, 
                Kbpsi_str.kbpsi_res_ion[ion2].data(), &num_prj, &zero, work_kbpsi, &num_prj);

        dgemm("N", "N", &num_prj, &num_orb, &num_prj, &mtwo, ddd, &num_prj, 
                Kbpsi_str.kbpsi_ion[ion2].data(), &num_prj, &one, work_kbpsi, &num_prj);


        for(st0 = 0; st0 < num_orb; st0++)
        {
            st1 = Kbpsi_str.orbital_index[ion2][st0];


            dgemv("N", &ct.max_nlpoints, &num_prj, &one, prjptr, &ct.max_nlpoints, 
                    &work_kbpsi[st0 * num_prj], &ione, &zero, prj_sum, &ione);


            /*
             *  project the prj_sum to the orbital |phi>  and stored in work 
             */



            qnm_beta_betapsi(&states1[st1], ion, prj_sum);


        }

        prjptr += pct.prj_per_ion[ion] * ct.max_nlpoints;       

    }                           /* end for ion */

    delete [] qnm_weight;
    delete [] prj_sum;
    delete [] work_kbpsi;
    

}
