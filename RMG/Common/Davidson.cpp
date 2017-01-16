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

#include <complex>
#include <omp.h>
#include <cmath>
#include <float.h>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "Mgrid.h"
#include "RmgException.h"
#include "Subdiag.h"
#include "Solvers.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "RmgParallelFft.h"
#include "TradeImages.h"
#include "packfuncs.h"

#include "transition.h"
#include "blas.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif



// Davidson diagonalization solver

template void Davidson<double>(Kpoint<double> *, double *, int &);
template void Davidson<std::complex<double> >(Kpoint<std::complex<double>> *, double *, int &);

#define DAVIDSON_DEBUG 0

static double occupied_tol = 0.01;

template <typename OrbitalType>
void Davidson (Kpoint<OrbitalType> *kptr, double *vtot, int &notconv)
{
    RmgTimer RT0("6-Davidson"), *RT1;

    OrbitalType alpha(1.0);
    OrbitalType beta(0.0);
    OrbitalType *newsint;
    if(ct.scf_steps > 0) {
        occupied_tol = std::min(occupied_tol, 0.1*ct.scf_accuracy / std::max(1.0, (double)ct.nel));
    }
    // Need this since the eigensolver may become unstable for very small residuals
    occupied_tol = std::max(occupied_tol, 1.0e-13);

    double unoccupied_tol = std::max(ct.unoccupied_tol_factor*occupied_tol, 1.0e-5 );

    //if(pct.gridpe == 0 && DAVIDSON_DEBUG)printf("OCCUPIED TOLERANCE = %20.12e\n",occupied_tol);

    int pbasis = kptr->pbasis;
    int nstates = kptr->nstates;
    notconv = nstates;
    int nbase = nstates;
    int max_steps = 10;
    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a;
    if(typeid(OrbitalType) == typeid(std::complex<double>)) {
         trans_a = trans_c;
    }
    else {
        trans_a = trans_t;
    }
    double vel = kptr->L->get_omega() / 
                 ((double)(kptr->G->get_NX_GRID(1) * kptr->G->get_NY_GRID(1) * kptr->G->get_NZ_GRID(1)));
    OrbitalType alphavel(vel);

    double avg_potential = 0.0;
    for(int idx = 0;idx < pbasis;idx++) avg_potential += vtot[idx];
    avg_potential = avg_potential / (double)pbasis;

    // For MPI routines
    int factor = 2;
    if(ct.is_gamma) factor = 1;
    OrbitalType *NULLptr = NULL;

    double *eigs = new double[ct.max_states];
    double *eigsw = new double[ct.max_states];
    bool *converged = new bool[ct.max_states]();


#if GPU_ENABLED
    OrbitalType *h_psi = (OrbitalType *)GpuMallocHost(pbasis * ct.max_states * sizeof(OrbitalType));
    OrbitalType *hr = (OrbitalType *)GpuMallocHost(ct.max_states * ct.max_states * sizeof(OrbitalType));
    OrbitalType *sr = (OrbitalType *)GpuMallocHost(ct.max_states * ct.max_states * sizeof(OrbitalType));
    OrbitalType *vr = (OrbitalType *)GpuMallocHost(ct.max_states * ct.max_states * sizeof(OrbitalType));
#else
    OrbitalType *h_psi = new OrbitalType[pbasis * ct.max_states];
    OrbitalType *hr = new OrbitalType[ct.max_states * ct.max_states]();
    OrbitalType *sr = new OrbitalType[ct.max_states * ct.max_states]();
    OrbitalType *vr = new OrbitalType[ct.max_states * ct.max_states]();
#endif

    for(int idx = 0;idx < nstates;idx++) vr[idx*ct.max_states + idx] = OrbitalType(1.0);

    // short version
    OrbitalType *psi = kptr->orbital_storage;

    // Apply Hamiltonian to current set of eigenvectors. At the current time
    // kptr->ns is also computed in ApplyHamiltonianBlock by AppNls and stored in kptr->ns
    RT1 = new RmgTimer("6-Davidson: Betaxpsi");
    Betaxpsi (kptr, 0, kptr->nstates, kptr->newsint_local, kptr->nl_weight);
    delete RT1;
    RT1 = new RmgTimer("6-Davidson: apply hamiltonian");
    double fd_diag = ApplyHamiltonianBlock (kptr, 0, nstates, h_psi, vtot); 
    delete RT1;
    OrbitalType *s_psi = kptr->ns;

    // Copy current eigs into compact array
    for(int st1 = 0;st1 < nstates;st1++) eigs[st1] = kptr->Kstates[st1].eig[0];

    // Compute A matrix
    RT1 = new RmgTimer("6-Davidson: matrix setup/reduce");
    RmgGemm(trans_a, trans_n, nbase, nbase, pbasis, alphavel, psi, pbasis, h_psi, pbasis, beta, hr, ct.max_states, 
            NULLptr, NULLptr, NULLptr, false, false, false, false);

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce it
    MPI_Request MPI_reqAij;
    MPI_Iallreduce(MPI_IN_PLACE, (double *)hr, nbase * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqAij);
#else
    MPI_Allreduce(MPI_IN_PLACE, (double *)hr, nbase * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
#endif

    // Compute S matrix
    RmgGemm (trans_a, trans_n, nbase, nbase, pbasis, alphavel, psi, pbasis, s_psi, pbasis, beta, sr, ct.max_states, 
             NULLptr, NULLptr, NULLptr, false, false, false, false);

#if HAVE_ASYNC_ALLREDUCE
    // Wait for Aij request to finish
    MPI_Wait(&MPI_reqAij, MPI_STATUS_IGNORE);
#endif

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce Sij request
    MPI_Request MPI_reqSij;
    MPI_Iallreduce(MPI_IN_PLACE, (double *)sr, nbase * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqSij);
#else
    MPI_Allreduce(MPI_IN_PLACE, (double *)sr, nbase * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
#endif

#if HAVE_ASYNC_ALLREDUCE
    // Wait for S request to finish and when done store copy in Sij
    MPI_Wait(&MPI_reqSij, MPI_STATUS_IGNORE);
#endif
    delete RT1;

    GeneralDiag(hr, sr, eigs, vr, nstates, nstates, ct.max_states, ct.subdiag_driver);
    for(int st=0;st < nstates;st++)kptr->Kstates[st].feig[0] = eigs[st];
    for(int st=0;st < nstates;st++)eigsw[st] = eigs[st];
    for(int st=0;st < nstates;st++)eigsw[st+nbase] = eigs[st];

    for(int steps = 0;steps < max_steps;steps++) {

        // Reorder eigenvectors
        int np = 0;
        for(int st = 0;st < nstates;st++) {

            if(!converged[st]) {

                if(np != st) {
                    for(int idx=0;idx < ct.max_states;idx++) vr[idx + np*ct.max_states] = vr[idx + st*ct.max_states];
                }
                eigsw[nbase + np] = eigs[st];
                np++;                

            }

        }

        // expand the basis set with the residuals ( H - e*S )|psi>
        RT1 = new RmgTimer("6-Davidson: generate residuals");
        RmgGemm(trans_n, trans_n, pbasis, notconv, nbase, alpha, s_psi, pbasis, vr, ct.max_states, beta, &psi[nbase*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, false, false, false);

#pragma omp parallel for
        for(int st1=0;st1 < notconv;st1++) {
            for(int idx=0;idx < pbasis;idx++) psi[(st1 + nbase)*pbasis + idx] = -eigsw[nbase + st1] * psi[(st1 + nbase)*pbasis + idx];
        }

        RmgGemm(trans_n, trans_n, pbasis, notconv, nbase, alpha, h_psi, pbasis, vr, ct.max_states, alpha, &psi[nbase*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, false, false, false);
        delete RT1;

        // Apply preconditioner
        RT1 = new RmgTimer("6-Davidson: precondition");
        DavPreconditioner (kptr, &psi[nbase*pbasis], fd_diag, &eigsw[nbase], vtot, notconv, avg_potential);
        delete RT1;

        // Normalize correction vectors. Not an exact normalization for norm conserving pseudopotentials
        // but that is OK. The goal is to get the magnitudes of all of the vectors being passed to the
        // diagonalizer roughly equal to improve stability.
        RT1 = new RmgTimer("6-Davidson: normalization");
        double *norms = new double[notconv]();
#pragma omp parallel for
        for(int st1=0;st1 < notconv;st1++) {
            for(int idx=0;idx < pbasis;idx++) norms[st1] += vel * std::norm(psi[(st1 + nbase)*pbasis + idx]);
        }

        MPI_Allreduce(MPI_IN_PLACE, (double *)norms, notconv, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

#pragma omp parallel for
        for(int st1=0;st1 < notconv;st1++) {
             norms[st1] = 1.0 / sqrt(norms[st1]);
             for(int idx=0;idx < pbasis;idx++) psi[(st1 + nbase)*pbasis + idx] *= norms[st1];
        }
        delete [] norms;
        delete RT1;


        // Apply Hamiltonian to the new vectors
        RT1 = new RmgTimer("6-Davidson: Betaxpsi");
        newsint = kptr->newsint_local + nbase*pct.num_nonloc_ions*ct.max_nl;
        Betaxpsi (kptr, nbase, notconv, newsint, kptr->nl_weight);
        delete RT1;

        RT1 = new RmgTimer("6-Davidson: apply hamiltonian");
        ApplyHamiltonianBlock (kptr, nbase, notconv, h_psi, vtot);
        delete RT1;


        // Update the reduced Hamiltonian and S matrices
        RT1 = new RmgTimer("6-Davidson: matrix setup/reduce");
        RmgGemm(trans_a, trans_n, nbase+notconv, notconv, pbasis, alphavel, psi, pbasis, &h_psi[nbase*pbasis], pbasis, beta, &hr[nbase*ct.max_states], ct.max_states, 
                NULLptr, NULLptr, NULLptr, false, false, false, false);

#if HAVE_ASYNC_ALLREDUCE
        // Asynchronously reduce it
        MPI_Request MPI_reqAij;
        MPI_Iallreduce(MPI_IN_PLACE, (double *)&hr[nbase*ct.max_states], notconv * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqAij);
#else
        MPI_Allreduce(MPI_IN_PLACE, (double *)&hr[nbase*ct.max_states], notconv * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
#endif

        RmgGemm(trans_a, trans_n, nbase+notconv, notconv, pbasis, alphavel, psi, pbasis, &s_psi[nbase*pbasis], pbasis, beta, &sr[nbase*ct.max_states], ct.max_states, 
                 NULLptr, NULLptr, NULLptr, false, false, false, false);

#if HAVE_ASYNC_ALLREDUCE
        // Wait for Aij request to finish
        MPI_Wait(&MPI_reqAij, MPI_STATUS_IGNORE);
#endif

#if HAVE_ASYNC_ALLREDUCE
        // Asynchronously reduce Sij request
        MPI_Request MPI_reqSij;
        MPI_Iallreduce(MPI_IN_PLACE, (double *)&sr[nbase*ct.max_states], notconv * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqSij);
#else
        MPI_Allreduce(MPI_IN_PLACE, (double *)&sr[nbase*ct.max_states], notconv * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
#endif

#if HAVE_ASYNC_ALLREDUCE
        // Wait for S request to finish
        MPI_Wait(&MPI_reqSij, MPI_STATUS_IGNORE);
#endif
        delete RT1;

        nbase = nbase + notconv;
        for(int i=0;i < nbase;i++) {
            for(int j=i+1;j < nbase;j++) {
                hr[j + i*ct.max_states] = hr[i + j*ct.max_states];
                sr[j + i*ct.max_states] = sr[i + j*ct.max_states];
            }
        }

        RT1 = new RmgTimer("6-Davidson: diagonalization");
        int info = GeneralDiag(hr, sr, eigsw, vr, nbase, nstates, ct.max_states, ct.subdiag_driver);
        delete RT1;
        if(info) {
            throw RmgFatalException() << "Diagonalization failed in Davidson, terminating." << " in " << __FILE__ << " at line " << __LINE__ << "\n";
        }

        // Check convergence
        int tnotconv = nstates;
        double max_tol = 0.0;
        double min_tol = DBL_MAX;
        int max_tol_state = 0;
        int min_tol_state = 0;
        double avg_occ_tol = 0.0;
        double avg_unocc_tol = 0.0;
        int occ_states = 0;
        int unocc_states = 0;
        double *tolerances = new double[nstates];

        for(int st=0;st < nstates;st++) {
            double occ = kptr->Kstates[st].occupation[0];
            if(ct.spin_flag) occ+= kptr->Kstates[st].occupation[1];
            tolerances[st] = fabs(eigs[st] - eigsw[st]);
            if((tolerances[st] > max_tol) && (occ > 0.002)) {
                max_tol = tolerances[st];
                max_tol_state = st;
                avg_occ_tol += tolerances[st];
                occ_states++;
            }
            if((tolerances[st] < min_tol) && (occ <= 0.002)) {
                min_tol = tolerances[st];
                min_tol_state = st;
                avg_unocc_tol += tolerances[st];
                unocc_states++;
            }
        }

        avg_occ_tol = avg_occ_tol / (double)occ_states;
        avg_unocc_tol = avg_unocc_tol / (double)unocc_states;

        for(int st=0;st < nstates;st++) {
            //printf("EIGS = %20.12f  %20.12f\n",eigs[st], eigsw[st]);
            double occ = kptr->Kstates[st].occupation[0];
            if(ct.spin_flag) occ+= kptr->Kstates[st].occupation[1];
            if(occ > 0.002) {
                converged[st] = (tolerances[st] < occupied_tol);
            }
            else {
                converged[st] = (tolerances[st] < unoccupied_tol);
            }
            if(converged[st]) tnotconv--;
            if((pct.gridpe==0) && (info != 0)) printf("STATE %d e0=%20.12f e1=%20.12f  TOLERANCE = %20.12f\n",st,eigs[st],eigsw[st],fabs(eigs[st] - eigsw[st]));
            //if((pct.gridpe==0)) printf("STATE %d e0=%20.12f e1=%20.12f  TOLERANCE = %20.12f\n",st,eigs[st],eigsw[st],fabs(eigs[st] - eigsw[st]));
        }
        delete [] tolerances;

        notconv = tnotconv;
        if(pct.gridpe==0 && DAVIDSON_DEBUG) rmg_printf("Davidson: notconv = %d  nbase=%d  occupied_tol=%7.3e\n", notconv, nbase, occupied_tol);
        if(pct.gridpe==0 && DAVIDSON_DEBUG) printf("MIN_TOLERANCE=%20.12e for STATE %d\n", min_tol, min_tol_state);
        if(pct.gridpe==0 && DAVIDSON_DEBUG) rmg_printf("MAX_TOLERANCE=%20.12e for STATE %d\n", max_tol, max_tol_state);
        if(pct.gridpe==0 && DAVIDSON_DEBUG) printf("AVG_OCC_TOLERANCE=%20.12e\n", avg_occ_tol);
        if(pct.gridpe==0 && DAVIDSON_DEBUG) printf("AVG_UNOCC_TOLERANCE=%20.12e\n", avg_unocc_tol);

        // Copy updated eigenvalues back
        for(int st=0;st < nstates;st++) eigs[st] = eigsw[st];

        // Check if we converged to the desired tolerance and if so return. If we
        // have exceeded the maximum number of iterations then we need to do something else.
        // If the expanded basis is getting too large then we need to rotate the orbitals
        // and start the davidson iteration again.
        if(((steps == (max_steps-1)) || ((nbase+notconv) > ct.max_states) || (notconv == 0))) {

            // Rotate orbitals
            RT1 = new RmgTimer("6-Davidson: rotate orbitals");
#if GPU_ENABLED
            OrbitalType *npsi = (OrbitalType *)GpuMallocHost(nstates*pbasis*sizeof(OrbitalType));
#else
            OrbitalType *npsi = new OrbitalType[nstates*pbasis];
#endif
            RmgGemm(trans_n, trans_n, pbasis, nstates, nbase, alpha, psi, pbasis, vr, ct.max_states, beta, npsi, pbasis, 
                NULLptr, NULLptr, NULLptr, false, false, false, false);
            for(int idx=0;idx < nstates*pbasis;idx++)psi[idx] = npsi[idx];
#if GPU_ENABLED
            GpuFreeHost(npsi);
#else
            delete [] npsi;
#endif
            delete RT1;

            if(notconv == 0) {
                rmg_printf("Davidson converged in %d steps\n", steps+1);
                break;  // done
            }

            if(steps == (max_steps-1)) {
                rmg_printf("Davidson incomplete convergence steps = %d\n", steps + 1);
                // Incomplete convergence, what should we do here?
                //throw RmgFatalException() << "Davidson failed to converge, terminating." << " in " << __FILE__ << " at line " << __LINE__ << "\n";
                break;

            }

            // refresh s_psi and h_psi
            RT1 = new RmgTimer("6-Davidson: refresh h_psi and s_psi");
            RmgGemm(trans_n, trans_n, pbasis, nstates, nbase, alpha, s_psi, pbasis, vr, ct.max_states, beta, &psi[nstates*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, false, false, false);
            for(int idx=0;idx < nstates*pbasis;idx++)s_psi[idx] = psi[nstates*pbasis + idx];

            RmgGemm(trans_n, trans_n, pbasis, nstates, nbase, alpha, h_psi, pbasis, vr, ct.max_states, beta, &psi[nstates*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, false, false, false);
            for(int idx=0;idx < nstates*pbasis;idx++)h_psi[idx] = psi[nstates*pbasis + idx];
            delete RT1;


            // Reset hr,sr,vr
            RT1 = new RmgTimer("6-Davidson: reset hr,sr,vr");
            nbase = nstates;
            for(int ix=0;ix < ct.max_states*ct.max_states;ix++) hr[ix] = OrbitalType(0.0);
            for(int ix=0;ix < ct.max_states*ct.max_states;ix++) sr[ix] = OrbitalType(0.0);
            for(int ix=0;ix < ct.max_states*ct.max_states;ix++) vr[ix] = OrbitalType(0.0);
            for(int st=0;st < nbase;st++) {
                hr[st + st*ct.max_states] = eigs[st];
                sr[st + st*ct.max_states] = OrbitalType(1.0);
                vr[st + st*ct.max_states] = OrbitalType(1.0);
            }
            delete RT1;
            
        }

    }

    // Copy eigs from compact array back into state structure
    for(int st = 0;st < nstates;st++) kptr->Kstates[st].eig[0] = eigs[st];



#if GPU_ENABLED
    GpuFreeHost(vr);
    GpuFreeHost(sr);
    GpuFreeHost(hr);
    GpuFreeHost(h_psi);
#else
    delete [] vr;
    delete [] sr;
    delete [] hr;
    delete [] h_psi;
#endif

    delete [] converged;
    delete [] eigsw;
    delete [] eigs;

    RT1 = new RmgTimer("6-Davidson: Betaxpsi");
    Betaxpsi (kptr, 0, kptr->nstates, kptr->newsint_local, kptr->nl_weight);
    delete RT1;
    
}

