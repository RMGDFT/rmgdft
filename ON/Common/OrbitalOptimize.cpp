/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>




#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"

#include "prototypes_on.h"
#include "init_var.h"
#include "transition.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"

#include "BaseThread.h"
#include "rmgthreads.h"
#include "RmgThread.h"
#include "PulayMixing.h"

extern std::vector<ORBITAL_PAIR> OrbitalPairs;

/* Flag for projector mixing */
static int mix_steps;

static void get_qnm_res(double *work_theta);
static void get_nonortho_res(STATE *, double *, STATE *);
static void get_dnmpsi(STATE *);

void OrbitalOptimize(STATE * states, STATE * states1, double *vxc, double *vh,
        double *vnuc, double *rho, double *rhoc, double * vxc_old, double * vh_old)
{
    int ione = 1;
    int order = ct.kohn_sham_fd_order;
    double hxgrid = Rmg_G->get_hxgrid(1);
    double hygrid = Rmg_G->get_hygrid(1);
    double hzgrid = Rmg_G->get_hzgrid(1);
    RmgTimer *RT = new RmgTimer("3-OrbitalOptimize");
    if(ct.LocalizedOrbitalLayout == LO_distribute)
    {
        int item = (ct.max_orbit_nx+order) *(ct.max_orbit_ny+order) *(ct.max_orbit_nz+order);
        FiniteDiff FD(&Rmg_L);



        RmgTimer *RT1a = new RmgTimer("3-OrbitalOptimize: distribute");
        DistributeToGlobal(vtot_c, vtot_global);
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
        MPI_Barrier(pct.grid_comm);
        delete(RT3);

        RmgTimer *RT4 = new RmgTimer("3-OrbitalOptimize: qnm");
        get_qnm_res(theta);
        MPI_Barrier(pct.grid_comm);
        delete(RT4);
        /* end shuchun wang */


        /* Loop over states istate to compute the whole matrix Hij 
           and all the vectors H|psi> */
        /* calculate the H |phi> on this processor and stored in states1[].psiR[] */

        RmgTimer *RTa = new RmgTimer("3-OrbitalOptimize: Hpsi");
        int st1;
        double *orbital_border;
        double *orbit_tem;
#pragma omp parallel private(st1,orbital_border,orbit_tem)
        {
            orbital_border = new double[2*item];
            orbit_tem = new double[2*item];
#pragma omp for schedule(static,1) nowait
            for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
            {
                STATE *sp = &states[st1];
                STATE *sp1 = &states1[st1];
                int ixx = states[st1].orbit_nx;
                int iyy = states[st1].orbit_ny;
                int izz = states[st1].orbit_nz;

                /* Generate 2*V*psi and store it  in orbit_tem */
                GenVxPsi(sp->psiR, st1, orbit_tem, vtot_global, states);

                double t1 = -1.0;
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

                if(ct.laplacian_offdiag || ct.laplacian_autocoeff)
                    FiniteDiffLap (orbital_border, orbit_tem, ixx, iyy, izz, LC);
                else
                    FD.app8_del2 (orbital_border, orbit_tem, ixx, iyy, izz, hxgrid, hygrid, hzgrid);

                t1 = 1.0;
                daxpy(&sp->size, &t1, orbit_tem, &ione, sp1->psiR, &ione);

            }                           /* end for st1 = .. */
            delete [] orbit_tem;
            delete [] orbital_border;
        }

        delete(RTa);
        RmgTimer *RT5 = new RmgTimer("3-OrbitalOptimize: dnm");
        get_dnmpsi(states1);
        delete(RT5);

    }


    /*
     *    precond(states1[ct.state_begin].psiR);
     *    t1 = 0.1; 
     *    daxpy(&pct.psi_size, &t1, states1[ct.state_begin].psiR, &ione, states[ct.state_begin].psiR, &ione); 
     */


    /*  SD, Pulay or KAIN method for Linear and Nonlinear equations */

    RmgTimer *RT6a = new RmgTimer("3-OrbitalOptimize: scale");

    delete(RT6a);
    if (ct.restart_mix == 1 || ct.move_centers_at_this_step == 1)
    {
        mix_steps = 0;
        Pulay_orbital->Refresh();
        if (pct.gridpe == 0)
            printf("\n restart the orbital mixing at step %d \n", ct.scf_steps);
    }


    RmgTimer *RT6 = new RmgTimer("3-OrbitalOptimize: mixing+precond");

    for (int st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        int ixx = states[st1].orbit_nx;
        int iyy = states[st1].orbit_ny;
        int izz = states[st1].orbit_nz;

        ZeroBoundary(states1[st1].psiR, ixx, iyy, izz);
        double residual = 0.0;
        for(int i = 0; i < ixx * iyy *izz; i++) residual += states1[st1].psiR[i] * states1[st1].psiR[i];
        //   printf("\n state residual %d %e", st1, residual);  

    }
    double gamma = -0.5;
    switch (ct.orbital_mixing_method)
    {
        case 0:
            Precond(states1[ct.state_begin].psiR);
            daxpy(&pct.psi_size, &gamma, states1[ct.state_begin].psiR, &ione, states[ct.state_begin].psiR, &ione);
            break;
        case 1:
            Pulay_orbital->Mixing(states[ct.state_begin].psiR, states1[ct.state_begin].psiR);
            break;
        case 2:
            Kain(mix_steps, pct.psi_size, states[ct.state_begin].psiR,
                    states1[ct.state_begin].psiR, ct.orbital_pulay_order);
            break;
        case 3:
            PulayWeighted(mix_steps, pct.psi_size, states[ct.state_begin].psiR,
                    states1[ct.state_begin].psiR, ct.orbital_pulay_order, 100, 0.5, 1);
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

    for (int st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        int ixx = states[st1].orbit_nx;
        int iyy = states[st1].orbit_ny;
        int izz = states[st1].orbit_nz;

        ZeroBoundary(states[st1].psiR, ixx, iyy, izz);
    }





#if     DEBUG
    print_status(states, vh, vxc, vnuc, vh_old, "before leaving OrbitalOptimize.c  ");
    print_state_sum(states1);
#endif
    delete(RT);

} 

/*-----------------------------------------------------------------*/

static void get_nonortho_res(STATE * states, double *work_theta, STATE * states1)
{
    int st1;

    //BaseThread *T = BaseThread::getBaseThread(0);

    //int active_threads = ct.MG_THREADS_PER_NODE;

    // theta * phi will be stored in states1.psiR, first set to 0.0
#pragma omp parallel private(st1)
    {
#pragma omp for schedule(static,1) nowait
        for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
            for (int idx = 0; idx < states[st1].size; idx++)
                states1[st1].psiR[idx] = 0.0;
    }

    if(OrbitalPairs.size() == 0 ) return;
    for(unsigned int ipair = 0; ipair < OrbitalPairs.size(); ipair++)
    {
        ORBITAL_PAIR onepair = OrbitalPairs[ipair];
        int st1 = onepair.orbital1;
        int st2 = onepair.orbital2;
        int st11 = st1-ct.state_begin;

        double temp = work_theta[st11 * ct.num_states + st2];
        ThetaPhi(st1, st2, temp, states[st2].psiR, states1[st1].psiR, 0, states, &onepair);
    }


    //  if( (int)OrbitalPairs.size() < active_threads) active_threads = (int)OrbitalPairs.size();

    //  int pair_start[active_threads], pair_end[active_threads];

    //  DistributeTasks(active_threads, (int)OrbitalPairs.size(), pair_start, pair_end);

    //  SCF_THREAD_CONTROL thread_control;

    //  for(int ist = 0;ist < active_threads;ist++) {



    //      thread_control.job = HYBRID_THETA_PHI;

    //      thread_control.nv = (void *)work_theta;

    //      thread_control.basetag = pair_start[ist];
    //      thread_control.extratag1 = active_threads;
    //      thread_control.extratag2 = pair_end[ist];

    //      QueueThreadTask(ist, thread_control);

    //      //DotProductOrbitOrbit(&states1[st1], &states[st2], &states[st1],  H, S, onepair);
    //  }

    //  // Thread tasks are set up so run them
    //  T->run_thread_tasks(active_threads);

}


void get_qnm_res(double *work_theta)
{

    int idx1;
    int num_prj, num_orb, tot_orb;
    int max_orb;
    double one = 1.0, zero = 0.0, *work_mat;

    max_orb = 0;

    for (int ion1 = 0; ion1 < pct.n_ion_center; ion1++)
    {
        tot_orb = Kbpsi_str.orbital_index[ion1].size();
        max_orb = std::max(max_orb, tot_orb);
    }

    work_mat = new double[(ct.state_end-ct.state_begin) *max_orb];


    for (int ion1 = 0; ion1 < pct.n_ion_center; ion1++)
    {

        num_prj = pct.prj_per_ion[pct.ionidx[ion1]];
        if (num_prj == 0) continue;
        num_orb = Kbpsi_str.num_orbital_thision[ion1]; 
        tot_orb = Kbpsi_str.orbital_index[ion1].size();

#pragma omp parallel private(idx1)
        {
#pragma omp for schedule(dynamic) nowait
            for(idx1 = 0; idx1 < num_orb; idx1++)
            {
                int st1 = Kbpsi_str.orbital_index[ion1][idx1];
                int st11 = st1 - ct.state_begin;

                for(int idx2 = 0; idx2 < tot_orb; idx2++)
                {
                    int st2 = Kbpsi_str.orbital_index[ion1][idx2];
                    work_mat[idx1 * tot_orb + idx2] = work_theta[st11 *ct.num_states + st2];
                }
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
    double *prjptr[pct.n_ion_center];
    double *prj_sum;

    double one=1.0, zero=0.0, mtwo=-2.0, *work_kbpsi; 

    /*
     *  dnm_weght[i] = sum_j dnm[i, j] *<beta_j|phi> 
     * prj_sum(r) = sum_i nm_weight[i] *beta_i(r)  
     */

    int num_orb_max = 1;
    for (int ion2 = 0; ion2 < pct.n_ion_center; ion2++)
    {
        num_orb_max = std::max(num_orb_max,Kbpsi_str.num_orbital_thision[ion2]); 
    }

    prj_sum = new double[num_orb_max * ct.max_nlpoints];
    work_kbpsi = new double[ct.max_nl * (ct.state_end-ct.state_begin)];


    size_t prj_ion_address = 0;
    for (int ion2 = 0; ion2 < pct.n_ion_center; ion2++)
    {
        int ion = pct.ionidx[ion2];
        prjptr[ion2] = &projectors[prj_ion_address];
        prj_ion_address += pct.prj_per_ion[ion] * ct.max_nlpoints;       
    }

    for (int ion2 = 0; ion2 < pct.n_ion_center; ion2++)
    {
        int ion = pct.ionidx[ion2];
        int num_prj = pct.prj_per_ion[ion];
        if (num_prj == 0) continue;
        double *qqq = pct.qqq[ion];
        double * ddd = pct.dnmI[ion];
        int num_orb = Kbpsi_str.num_orbital_thision[ion2]; 

        dgemm("N", "N", &num_prj, &num_orb, &num_prj, &one, qqq, &num_prj, 
                Kbpsi_str.kbpsi_res_ion[ion2].data(), &num_prj, &zero, work_kbpsi, &num_prj);

        dgemm("N", "N", &num_prj, &num_orb, &num_prj, &mtwo, ddd, &num_prj, 
                Kbpsi_str.kbpsi_ion[ion2].data(), &num_prj, &one, work_kbpsi, &num_prj);

        dgemm("N", "N", &ct.max_nlpoints, &num_orb, &num_prj, &one, prjptr[ion2], &ct.max_nlpoints, 
                work_kbpsi, &num_prj, &zero, prj_sum, &ct.max_nlpoints);
        int st0;
#pragma omp parallel private(st0)
        {
            //#pragma omp for schedule(static,1) nowait
#pragma omp for schedule(dynamic,1) nowait
            for(st0 = 0; st0 < num_orb; st0++)
            {
                int st1 = Kbpsi_str.orbital_index[ion2][st0];


                /*
                 *  project the prj_sum to the orbital |phi>  and stored in work 
                 */


                qnm_beta_betapsi(&states1[st1], ion, &prj_sum[ct.max_nlpoints * st0]);


            }

        }                           /* end for ion */

    }
    delete [] prj_sum;
    delete [] work_kbpsi;


}




