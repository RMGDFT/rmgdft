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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include <complex>
#include "RmgParallelFft.h"
#include "GlobalSums.h"
#include "Prolong.h"
#include "rmgthreads.h"
#include "RmgThread.h"
#include "Symmetry.h"
#include "Voronoi.h"



template void GetNewRho<double>(Kpoint<double> **, double *);
template void GetNewRho<std::complex<double> >(Kpoint<std::complex<double>> **, double *);

template void GetNewRhoPre<double>(Kpoint<double> **, double *);
template void GetNewRhoPre<std::complex<double> >(Kpoint<std::complex<double>> **, double *);

template void GetNewRhoPost<double>(Kpoint<double> **, double *);
template void GetNewRhoPost<std::complex<double> >(Kpoint<std::complex<double>> **, double *);

template void GetNewRhoOne<double>(double *, Prolong *, double *, double);
template void GetNewRhoOne<std::complex<double>>(std::complex<double> *, Prolong *, double *, double);


template <typename OrbitalType> void GetNewRho(Kpoint<OrbitalType> **Kpts, double *rho)
{
    int factor = ct.noncoll_factor * ct.noncoll_factor;
    int ratio = Rmg_G->default_FG_RATIO;
    int FP0_BASIS = Rmg_G->get_P0_BASIS(ratio);

    if(ct.fast_density)
        GetNewRhoPost(Kpts, rho);
    else
        GetNewRhoPre(Kpts, rho);
        
    if(!ct.norm_conserving_pp) {
        double *augrho = new double[FP0_BASIS*factor]();
        GetAugRho(Kpts, augrho);
        for(int idx = 0;idx < FP0_BASIS*factor;idx++) rho[idx] += augrho[idx];
        delete [] augrho;
    }


    if(Rmg_Symm) Rmg_Symm->symmetrize_grid_object(rho);
    if(ct.noncoll && Rmg_Symm)
        Rmg_Symm->symmetrize_grid_vector(&rho[FP0_BASIS]);


    /* Check total charge. */
    ct.tcharge = ZERO;
    for (int idx = 0; idx < FP0_BASIS; idx++)
        ct.tcharge += rho[idx];

    /* ct.tcharge = real_sum_all (ct.tcharge); */
    ct.tcharge = real_sum_all (ct.tcharge, pct.grid_comm);
    ct.tcharge = real_sum_all (ct.tcharge, pct.spin_comm);
    ct.tcharge = ct.tcharge * get_vel_f();

    /* Renormalize charge, there could be some discrepancy because of interpolation */
    double t1 = ct.nel / ct.tcharge;
    for(int i = 0;i < FP0_BASIS * factor;i++) rho[i] *= t1;

    /*Write out normalization constant if needed*/
    double difference = fabs(t1 - 1.0);
    if ((ct.verbose == 1) || (difference > 0.01))
        rmg_printf ("Charge normalization constant: %f\n", t1);

}

// Generates the new density by interpolating each orbital to the fine grid and then squaring
// and summing them. 
template <typename OrbitalType> void GetNewRhoPre(Kpoint<OrbitalType> **Kpts, double *rho)
{

    if(Verify ("freeze_occupied", true, Kpts[0]->ControlMap)) return;

    BaseThread *T = BaseThread::getBaseThread(0);
    int nstates = Kpts[0]->nstates;
    int ratio = Rmg_G->default_FG_RATIO;
    int FP0_BASIS = Rmg_G->get_P0_BASIS(ratio);
    Prolong P(ratio, 10, *Rmg_T);

    int factor = ct.noncoll_factor * ct.noncoll_factor;
    double *work = new double[FP0_BASIS * factor * ct.MG_THREADS_PER_NODE]();


    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode) active_threads--;
    if(active_threads < 1) active_threads = 1;

    int istop = nstates / active_threads;
    istop = istop * active_threads;

    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {

        /* Loop over states and accumulate charge */
        for(int st1=0;st1 < istop;st1+=active_threads)
        {

            SCF_THREAD_CONTROL thread_control;

            for(int ist = 0;ist < active_threads;ist++)
            {
                if(fabs(Kpts[kpt]->Kstates[st1+ist].occupation[0]) > 1.0e-10)
                {
		    thread_control.job = HYBRID_GET_RHO;
                    double scale = Kpts[kpt]->Kstates[st1+ist].occupation[0] * Kpts[kpt]->kp.kweight;
                    OrbitalType *psi = Kpts[kpt]->Kstates[st1+ist].psi;
                    thread_control.p1 = (void *)psi;
                    thread_control.p2 = (void *)&P;
                    thread_control.p3 = (void *)&work[ist*FP0_BASIS * factor];
                    thread_control.fd_diag = scale;
                    thread_control.basetag = st1 + ist;
                    QueueThreadTask(ist, thread_control);
                }
            }
            // Thread tasks are set up so wake them
            if(!ct.mpi_queue_mode) T->run_thread_tasks(active_threads);

        } 
        if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

        // Process any remaining states in serial fashion
        for(int st1 = istop;st1 < nstates;st1++) 
        {
            if(fabs(Kpts[kpt]->Kstates[st1].occupation[0]) > 1.0e-10)
            {
                double scale = Kpts[kpt]->Kstates[st1].occupation[0] * Kpts[kpt]->kp.kweight;
                OrbitalType *psi = Kpts[kpt]->Kstates[st1].psi;
                GetNewRhoOne(psi, &P, work, scale);
            }
        }

    }                           /*end for kpt */

    // Combine contributions from all of the threads
    for(int st2=1;st2 < active_threads;st2++)
    {
        for(int idx=0;idx < FP0_BASIS*factor;idx++)
        {
            work[idx] += work[st2*FP0_BASIS*factor+idx];
        } 
    }

    MPI_Allreduce(MPI_IN_PLACE, (double *)work, FP0_BASIS * factor, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
    if(ct.noncoll)
    {
        double rho_up, rho_down;
        for(int idx = 0; idx < FP0_BASIS; idx++)
        {
            rho_up = work[idx];
            rho_down = work[idx+3 * FP0_BASIS];
            work[idx] = rho_up + rho_down;
            work[idx+3*FP0_BASIS] = rho_up - rho_down;
        }
    }

    for(int idx = 0; idx < FP0_BASIS * factor; idx++) rho[idx] = work[idx];
    delete [] work;
}



template <typename OrbitalType> void GetNewRhoOne(OrbitalType *psi, Prolong *P, double *work, double scale)
{

    int ratio = Rmg_G->default_FG_RATIO;
    int FP0_BASIS = Rmg_G->get_P0_BASIS(ratio);
    int dimx = Rmg_G->get_PX0_GRID(ratio);
    int dimy = Rmg_G->get_PY0_GRID(ratio);
    int dimz = Rmg_G->get_PZ0_GRID(ratio);
    int half_dimx = Rmg_G->get_PX0_GRID(1);
    int half_dimy = Rmg_G->get_PY0_GRID(1);
    int half_dimz = Rmg_G->get_PZ0_GRID(1);

    std::complex<double> psiud;
    OrbitalType *psi_f = new OrbitalType[ct.noncoll_factor * FP0_BASIS]();

    P->prolong(psi_f, psi, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);
    if(ct.noncoll)
       P->prolong(&psi_f[FP0_BASIS], &psi[half_dimx*half_dimy*half_dimz], dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);

    for (int idx = 0; idx < FP0_BASIS; idx++)
    {
        work[idx] += scale * std::norm(psi_f[idx]);
        if(ct.noncoll)
        {
            psiud = 2.0 * psi_f[idx] * std::conj(psi_f[idx + FP0_BASIS]);
            work[idx + 1 * FP0_BASIS] += scale * std::real(psiud);
            work[idx + 2 * FP0_BASIS] += scale * std::imag(psiud);
            work[idx + 3 * FP0_BASIS] += scale * std::norm(psi_f[idx + FP0_BASIS]);
        }
    }                   /* end for */

    delete [] psi_f;

}


// Computes the sum of the squares of the orbitals on the coarse grid and then does a single
// interpolation of the resulting charge density to the fine grid. Faster but less accurate.

template <typename OrbitalType> void GetNewRhoPost(Kpoint<OrbitalType> **Kpts, double *rho)
{

    int pbasis = Kpts[0]->pbasis;
    int nstates = Kpts[0]->nstates;

    if(Verify ("freeze_occupied", true, Kpts[0]->ControlMap)) return;

    std::complex<double> psiud;
    int factor = ct.noncoll_factor * ct.noncoll_factor;
    double *work = new double[pbasis * factor];

    for(int idx = 0;idx < pbasis * factor;idx++)
        work[idx] = 0.0;

    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
#pragma omp parallel
        {
            double *tarr = new double[pbasis * factor]();
#pragma omp barrier

            /* Loop over states and accumulate charge */
#pragma omp for schedule(static, 1) nowait
            for (int istate = 0; istate < nstates; istate++)
            {

                double scale = Kpts[kpt]->Kstates[istate].occupation[0] * Kpts[kpt]->kp.kweight;

                OrbitalType *psi = Kpts[kpt]->Kstates[istate].psi;

                for (int idx = 0; idx < pbasis; idx++)
                {
                    tarr[idx] += scale * std::norm(psi[idx]);
                    if(ct.noncoll)
                    {
                        psiud = 2.0 * psi[idx] * std::conj(psi[idx + pbasis]);
                        tarr[idx + 1 * pbasis] += scale * std::real(psiud);
                        tarr[idx + 2 * pbasis] += scale * std::imag(psiud);
                        tarr[idx + 3 * pbasis] += scale * std::norm(psi[idx + pbasis]);
                    }
                }                   /* end for */

            }                       /*end for istate */
#pragma omp critical
            for(int idx = 0; idx < pbasis * factor; idx++) work[idx] += tarr[idx];
            delete [] tarr;
        }
    }                           /*end for kpt */

    MPI_Allreduce(MPI_IN_PLACE, (double *)work, pbasis * factor, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
    if(ct.noncoll)
    {
        double rho_up, rho_down;
        for(int idx = 0; idx < pbasis; idx++)
        {
            rho_up = work[idx];
            rho_down = work[idx+3 * pbasis];
            work[idx] = rho_up + rho_down;
            work[idx+3*pbasis] = rho_up - rho_down;
        }
    }


    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    /* Interpolate onto fine grid, result will be stored in rho*/
    for(int is = 0; is < factor; is++)
    {
        switch (ct.interp_flag)
        {
            case CUBIC_POLYNOMIAL_INTERPOLATION:
                pack_rho_ctof (&work[is*pbasis], &rho[is*FP0_BASIS]);
                break;
            case PROLONG_INTERPOLATION:
                mg_prolong_MAX10 (&rho[is*FP0_BASIS], &work[is*pbasis], get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);
                break;
            case FFT_INTERPOLATION:
                FftInterpolation (*Kpts[0]->G, &work[is*pbasis], &rho[is*FP0_BASIS], Rmg_G->default_FG_RATIO, ct.sqrt_interpolation);
                break;
            default:

                //Dprintf ("charge interpolation is set to %d", ct.interp_flag);
                rmg_error_handler (__FILE__, __LINE__, "ct.interp_flag is set to an invalid value.");


        }
    }

    delete [] work;
}
