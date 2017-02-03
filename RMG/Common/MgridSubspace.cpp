/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "BaseThread.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Subdiag.h"
#include "rmgthreads.h"
#include "packfuncs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "../Headers/prototypes.h"


// Solver that uses multigrid preconditioning and subspace rotations

template <typename OrbitalType> void MgridSubspace (Kpoint<OrbitalType> *kptr, double *);



template void MgridSubspace<double> (Kpoint<double> *, double *);
template void MgridSubspace<std::complex<double> > (Kpoint<std::complex<double>> *, double *);

template <typename OrbitalType> void MgridSubspace (Kpoint<OrbitalType> *kptr, double *vtot_psi)
{
    RmgTimer RT0("3-MgridSubspace"), *RT1;
    BaseThread *T = BaseThread::getBaseThread(0);

    double mean_occ_res = DBL_MAX;
    double mean_unocc_res = DBL_MAX;
    double max_occ_res = 0.0;
    double max_unocc_res = 0.0;
    double min_occ_res = DBL_MAX;
    double min_unocc_res = DBL_MAX;
    bool potential_acceleration = ((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0));
    int pbasis = kptr->pbasis;

    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode) active_threads--;
    if(active_threads < 1) active_threads = 1;

    for(int vcycle = 0;vcycle < ct.eig_parm.mucycles;vcycle++) {

        // Update betaxpsi        
        RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
        Betaxpsi (kptr, 0, kptr->nstates, kptr->newsint_local, kptr->nl_weight);
        delete(RT1);

        /* Update the wavefunctions */
        int istop = kptr->nstates / active_threads;
        istop = istop * active_threads;

        // Apply the non-local operators to a block of orbitals
        RT1 = new RmgTimer("3-MgridSubspace: AppNls");
        AppNls(kptr, kptr->newsint_local, kptr->Kstates[0].psi, kptr->nv, kptr->ns, kptr->Bns,
               0, std::min(ct.non_local_block_size, kptr->nstates));
        delete(RT1);
        int first_nls = 0;

        for(int st1=0;st1 < istop;st1+=active_threads) {
          SCF_THREAD_CONTROL thread_control;

          // Make sure the non-local operators are applied for the next block if needed
          int check = first_nls + active_threads;
          if(check > ct.non_local_block_size) {
              RT1 = new RmgTimer("3-MgridSubspace: AppNls");
              AppNls(kptr, kptr->newsint_local, kptr->Kstates[st1].psi, kptr->nv, &kptr->ns[st1 * pbasis], kptr->Bns,
                     st1, std::min(ct.non_local_block_size, kptr->nstates - st1));
              first_nls = 0;
              delete(RT1);
          }
        
          RT1 = new RmgTimer("3-MgridSubspace: Mg_eig");
          for(int ist = 0;ist < active_threads;ist++) {
              thread_control.job = HYBRID_EIG;
              thread_control.vtot = vtot_psi;
              thread_control.vcycle = vcycle;
              thread_control.sp = &kptr->Kstates[st1 + ist];
              thread_control.p3 = (void *)kptr;
              thread_control.nv = (void *)&kptr->nv[(first_nls + ist) * pbasis];
              thread_control.ns = (void *)&kptr->ns[(st1 + ist) * pbasis];  // ns is not blocked!
              thread_control.basetag = kptr->Kstates[st1 + ist].istate;
              QueueThreadTask(ist, thread_control);
          }

          // Thread tasks are set up so run them
          if(!ct.mpi_queue_mode) T->run_thread_tasks(active_threads);

          // Increment index into non-local block
          first_nls += active_threads;
          delete RT1;

        }

        if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

        // Process any remaining states in serial fashion
        RT1 = new RmgTimer("3-MgridSubspace: Mg_eig");
        for(int st1 = istop;st1 < kptr->nstates;st1++) {
            if(ct.is_gamma) {
                if(ct.rms > ct.preconditioner_thr)
                    MgEigState<double,float> ((Kpoint<double> *)kptr, (State<double> *)&kptr->Kstates[st1], vtot_psi, 
                                               (double *)&kptr->nv[first_nls * pbasis], (double *)&kptr->ns[st1 * pbasis], vcycle);
                else
                    MgEigState<double,double> ((Kpoint<double> *)kptr, (State<double> *)&kptr->Kstates[st1], vtot_psi, 
                                               (double *)&kptr->nv[first_nls * pbasis], (double *)&kptr->ns[st1 * pbasis], vcycle);
            }
            else {
                if(ct.rms > ct.preconditioner_thr)
                    MgEigState<std::complex<double>, std::complex<float> > ((Kpoint<std::complex<double>> *)kptr, (State<std::complex<double> > *)&kptr->Kstates[st1], vtot_psi, 
                                               (std::complex<double> *)&kptr->nv[first_nls * pbasis], (std::complex<double> *)&kptr->ns[st1 * pbasis], vcycle);
                else
                    MgEigState<std::complex<double>, std::complex<double> > ((Kpoint<std::complex<double>> *)kptr, (State<std::complex<double> > *)&kptr->Kstates[st1], vtot_psi, 
                                               (std::complex<double> *)&kptr->nv[first_nls * pbasis], (std::complex<double> *)&kptr->ns[st1 * pbasis], vcycle);
            }
            first_nls++;
        }
        delete(RT1);

    }

    if(Verify ("freeze_occupied", true, kptr->ControlMap)) {

        // Orbital residual measures (used for some types of calculations
        kptr->max_unocc_res_index = (int)(ct.gw_residual_fraction * (double)kptr->nstates);
        kptr->mean_occ_res = 0.0;
        kptr->min_occ_res = DBL_MAX;
        kptr->max_occ_res = 0.0;
        kptr->mean_unocc_res = 0.0;
        kptr->min_unocc_res = DBL_MAX;
        kptr->max_unocc_res = 0.0;
        kptr->highest_occupied = 0;
        for(int istate = 0;istate < kptr->nstates;istate++) {
            if(kptr->Kstates[istate].occupation[0] > 0.0) {
                kptr->mean_occ_res += kptr->Kstates[istate].res;
                mean_occ_res += kptr->Kstates[istate].res;
                if(kptr->Kstates[istate].res >  kptr->max_occ_res)  kptr->max_occ_res = kptr->Kstates[istate].res;
                if(kptr->Kstates[istate].res <  kptr->min_occ_res)  kptr->min_occ_res = kptr->Kstates[istate].res;
                if(kptr->Kstates[istate].res >  max_occ_res)  max_occ_res = kptr->Kstates[istate].res;
                if(kptr->Kstates[istate].res <  min_occ_res)  min_occ_res = kptr->Kstates[istate].res;
                kptr->highest_occupied = istate;
            }
            else {
                if(istate <= kptr->max_unocc_res_index) {
                    kptr->mean_unocc_res += kptr->Kstates[istate].res;
                    mean_unocc_res += kptr->Kstates[istate].res;
                    if(kptr->Kstates[istate].res >  kptr->max_unocc_res)  kptr->max_unocc_res = kptr->Kstates[istate].res;
                    if(kptr->Kstates[istate].res <  kptr->min_unocc_res)  kptr->min_unocc_res = kptr->Kstates[istate].res;
                    if(kptr->Kstates[istate].res >  max_unocc_res)  max_unocc_res = kptr->Kstates[istate].res;
                    if(kptr->Kstates[istate].res <  min_unocc_res)  min_unocc_res = kptr->Kstates[istate].res;
                }
            }
        }
        kptr->mean_occ_res = kptr->mean_occ_res / (double)(kptr->highest_occupied + 1);
        kptr->mean_unocc_res = kptr->mean_unocc_res / (double)(kptr->max_unocc_res_index -(kptr->highest_occupied + 1));
        mean_occ_res = mean_occ_res / (double)(ct.num_kpts*(kptr->highest_occupied + 1));
        mean_unocc_res = mean_unocc_res / (double)(ct.num_kpts*kptr->max_unocc_res_index -(kptr->highest_occupied + 1));

        rmg_printf("Mean/Min/Max unoccupied wavefunction residual for kpoint %d  =  %10.5e  %10.5e  %10.5e\n", kptr->kidx, kptr->mean_unocc_res, kptr->min_unocc_res, kptr->max_unocc_res);

    }


    /* wavefunctions have changed, projectors have to be recalculated
     * but if we are using potential acceleration and not well converged yet
     * it is counterproductive to do so */
    if(!potential_acceleration || (potential_acceleration && (ct.rms <  5.0e-6))) {
        RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
        Betaxpsi (kptr, 0, kptr->nstates, kptr->newsint_local, kptr->nl_weight);

        delete(RT1);
    }


    /* Now we orthognalize and optionally do subspace diagonalization
     * In the gamma point case, orthogonalization is not required when doing subspace diagonalization
     * For non-gamma point we have to do first orthogonalization and then, optionally subspace diagonalization
     * the reason is for non-gamma subdiag is not coded to solve generalized eigenvalue problem, it can
     * only solve the regular eigenvalue problem and that requires that wavefunctions are orthogonal to start with.*/

    bool diag_this_step = (ct.diag && ct.scf_steps % ct.diag == 0 && ct.scf_steps < ct.end_diag);

    /* do diagonalizations if requested, if not orthogonalize */
    if (diag_this_step) {

        RT1 = new RmgTimer("3-MgridSubspace: Diagonalization");
        Subdiag (kptr, vtot_psi, ct.subdiag_driver);
        delete(RT1);
        RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
        Betaxpsi (kptr, 0, kptr->nstates, kptr->newsint_local, kptr->nl_weight);
        delete(RT1);
        // Projectors are rotated along with orbitals in Subdiag so no need to recalculate
        // after diagonalizing.

    }
    else {

        RT1 = new RmgTimer("3-MgridSubspace: Orthogonalization");
        kptr->orthogonalize(kptr->orbital_storage);
        delete(RT1);

        // wavefunctions have changed, projectors have to be recalculated */
        RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
        Betaxpsi (kptr, 0, kptr->nstates, kptr->newsint_local, kptr->nl_weight);
        delete(RT1);

    }
        

    /* If sorting is requested then sort the states. */
    if (ct.sortflag) {
        kptr->sort_orbitals();
    }


}


