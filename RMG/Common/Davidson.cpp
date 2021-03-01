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



// Davidson diagonalization solver part of Kpoint class
template void Kpoint<double>::Davidson(double *vtot, double *vxc_psi, int &notconv);
template void Kpoint<std::complex<double>>::Davidson(double *vtot, double *vxc_psi, int &notconv);


#define DAVIDSON_DEBUG 0

static double occupied_tol = 0.01;

template <class KpointType> void Kpoint<KpointType>::Davidson(double *vtot, double *vxc_psi, int &notconv)
{
    RmgTimer RT0("6-Davidson"), *RT1;

    KpointType alpha(1.0);
    KpointType beta(0.0);
    KpointType *newsint;

    KpointType *weight = this->nl_weight;
#if HIP_ENABLED
    weight = this->nl_weight_gpu;
#endif

    double acheck = ct.scf_accuracy;
    if(ct.scf_steps == 0) acheck = 0.01;
    occupied_tol = 0.1*acheck / std::max(1.0, (double)ct.nel);
    if(ct.spinorbit || ct.noncoll) occupied_tol /= 8.0;

    occupied_tol = std::min(occupied_tol, 1.0e-4);
    // Need this since the eigensolver may become unstable for very small residuals
    occupied_tol = std::max(occupied_tol, 1.0e-13);
    double unoccupied_tol = std::max(ct.unoccupied_tol_factor*occupied_tol, 1.0e-4 );
    if(ct.spinorbit || ct.noncoll) unoccupied_tol /= 8.0;

    //if(pct.gridpe == 0 && DAVIDSON_DEBUG)printf("OCCUPIED TOLERANCE = %20.12e\n",occupied_tol);

    notconv = nstates;
    int nbase = nstates;
    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a;
    if(typeid(KpointType) == typeid(std::complex<double>)) {
         trans_a = trans_c;
    }
    else {
        trans_a = trans_t;
    }
    double vel = this->L->get_omega() / 
                 ((double)(this->G->get_NX_GRID(1) * this->G->get_NY_GRID(1) * this->G->get_NZ_GRID(1)));
    KpointType alphavel(vel);

    double avg_potential = 0.0;
    for(int idx = 0;idx < pbasis;idx++) avg_potential += vtot[idx];
    avg_potential = avg_potential / (double)pbasis;

    // For MPI routines
    int factor = 2;
    if(ct.is_gamma) factor = 1;

    double *eigs = new double[ct.max_states];
    double *eigsw = new double[2*ct.max_states];
    bool *converged = new bool[ct.max_states]();

#if CUDA_ENABLED || HIP_ENABLED
    KpointType *h_psi = (KpointType *)RmgMallocHost(pbasis_noncoll * ct.max_states * sizeof(KpointType));
    KpointType *hr = (KpointType *)RmgMallocHost(ct.max_states * ct.max_states * sizeof(KpointType));
    KpointType *sr = (KpointType *)RmgMallocHost(ct.max_states * ct.max_states * sizeof(KpointType));
    KpointType *vr = (KpointType *)RmgMallocHost(ct.max_states * ct.max_states * sizeof(KpointType));
#else
    KpointType *h_psi = new KpointType[pbasis_noncoll * ct.max_states];
    KpointType *hr = new KpointType[ct.max_states * ct.max_states]();
    KpointType *sr = new KpointType[ct.max_states * ct.max_states]();
    KpointType *vr = new KpointType[ct.max_states * ct.max_states]();
#endif

    for(int idx = 0;idx < nstates;idx++) vr[idx*ct.max_states + idx] = KpointType(1.0);

    // short version
    KpointType *psi = this->orbital_storage;

    // Apply Hamiltonian to current set of eigenvectors. At the current time
    // this->ns is also computed in ApplyHamiltonianBlock by AppNls and stored in this->ns
    RT1 = new RmgTimer("6-Davidson: Betaxpsi");
    this->BetaProjector->project(this, this->newsint_local, 0, nstates*ct.noncoll_factor, weight);
    delete RT1;

    if(ct.ldaU_mode != LDA_PLUS_U_NONE)
    {   
        RmgTimer RTL("6-Davidson: ldaUop x psi"); 
        LdaplusUxpsi(this, 0, this->nstates, this->orbitalsint_local);
    }

    RT1 = new RmgTimer("6-Davidson: apply hamiltonian");
    double fd_diag = ApplyHamiltonianBlock<KpointType> (this, 0, nstates, h_psi, vtot, vxc_psi); 
    delete RT1;
    KpointType *s_psi = this->ns;
    if(ct.norm_conserving_pp && ct.is_gamma) s_psi = this->orbital_storage;

    // Copy current eigs into compact array
    for(int st1 = 0;st1 < nstates;st1++) eigs[st1] = this->Kstates[st1].eig[0];

    // Compute A matrix
    RT1 = new RmgTimer("6-Davidson: matrix setup/reduce");
    RmgGemm(trans_a, trans_n, nbase, nbase, pbasis_noncoll, alphavel, psi, pbasis_noncoll, h_psi, pbasis_noncoll, beta, hr, ct.max_states);

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce it
    MPI_Request MPI_reqAij;
    if(ct.use_async_allreduce)
        MPI_Iallreduce(MPI_IN_PLACE, (double *)hr, nbase * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqAij);
    else
        BlockAllreduce((double *)hr, (size_t)nbase*(size_t)ct.max_states * (size_t)factor, pct.grid_comm);
#else
    BlockAllreduce((double *)hr, (size_t)nbase*(size_t)ct.max_states * (size_t)factor, pct.grid_comm);

#endif

    // Compute S matrix
    RmgGemm (trans_a, trans_n, nbase, nbase, pbasis_noncoll, alphavel, psi, pbasis_noncoll, s_psi, pbasis_noncoll, beta, sr, ct.max_states);

#if HAVE_ASYNC_ALLREDUCE
    // Wait for Aij request to finish
    if(ct.use_async_allreduce) MPI_Wait(&MPI_reqAij, MPI_STATUS_IGNORE);
#endif

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce Sij request
    MPI_Request MPI_reqSij;
    if(ct.use_async_allreduce)
        MPI_Iallreduce(MPI_IN_PLACE, (double *)sr, nbase * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqSij);
    else
        BlockAllreduce((double *)sr, (size_t)nbase*(size_t)ct.max_states * (size_t)factor, pct.grid_comm);
#else
    BlockAllreduce((double *)sr, (size_t)nbase*(size_t)ct.max_states * (size_t)factor, pct.grid_comm);
#endif

#if HAVE_ASYNC_ALLREDUCE
    // Wait for S request to finish and when done store copy in Sij
    if(ct.use_async_allreduce) MPI_Wait(&MPI_reqSij, MPI_STATUS_IGNORE);
#endif
    delete RT1;

    GeneralDiag(hr, sr, eigs, vr, nstates, nstates, ct.max_states, ct.subdiag_driver);
    for(int st=0;st < nstates;st++)this->Kstates[st].feig[0] = eigs[st];
    for(int st=0;st < nstates;st++)this->Kstates[st].eig[0] = eigs[st];
    for(int st=0;st < nstates;st++)eigsw[st] = eigs[st];
    for(int st=0;st < nstates;st++)eigsw[st+nbase] = eigs[st];

    for(int steps = 0;steps < ct.david_max_steps;steps++) {

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
        RmgGemm(trans_n, trans_n, pbasis_noncoll, notconv, nbase, alpha, s_psi, pbasis_noncoll, vr, ct.max_states, beta, &psi[nbase*pbasis_noncoll], pbasis_noncoll);

#pragma omp parallel for
        for(int st1=0;st1 < notconv;st1++) {
            for(int idx=0;idx < pbasis_noncoll;idx++) psi[(st1 + nbase)*pbasis_noncoll + idx] = -eigsw[nbase + st1] * psi[(st1 + nbase)*pbasis_noncoll + idx];
        }

        RmgGemm(trans_n, trans_n, pbasis_noncoll, notconv, nbase, alpha, h_psi, pbasis_noncoll, vr, ct.max_states, alpha, &psi[nbase*pbasis_noncoll], pbasis_noncoll);
        delete RT1;

        // Apply preconditioner
        RT1 = new RmgTimer("6-Davidson: precondition");
        DavPreconditioner (this, &psi[nbase*pbasis_noncoll], fd_diag, &eigsw[nbase], vtot, notconv, avg_potential);
        delete RT1;

        // Normalize correction vectors. Not an exact normalization for norm conserving pseudopotentials
        // but that is OK. The goal is to get the magnitudes of all of the vectors being passed to the
        // diagonalizer roughly equal to improve stability.
        RT1 = new RmgTimer("6-Davidson: normalization");
        double *norms = new double[notconv]();
#pragma omp parallel for
        for(int st1=0;st1 < notconv;st1++) {
            for(int idx=0;idx < pbasis_noncoll;idx++) norms[st1] += vel * std::norm(psi[(st1 + nbase)*pbasis_noncoll + idx]);
        }

        MPI_Allreduce(MPI_IN_PLACE, (double *)norms, notconv, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

#pragma omp parallel for
        for(int st1=0;st1 < notconv;st1++) {
             norms[st1] = 1.0 / sqrt(norms[st1]);
             for(int idx=0;idx < pbasis_noncoll;idx++) psi[(st1 + nbase)*pbasis_noncoll + idx] *= norms[st1];
        }
        delete [] norms;
        delete RT1;


        // Apply Hamiltonian to the new vectors
        RT1 = new RmgTimer("6-Davidson: Betaxpsi");
        newsint = this->newsint_local + nbase * this->BetaProjector->get_num_nonloc_ions() * ct.max_nl * ct.noncoll_factor;
        this->BetaProjector->project(this, newsint, nbase*ct.noncoll_factor, notconv*ct.noncoll_factor, weight);
        delete RT1;

        if(ct.ldaU_mode != LDA_PLUS_U_NONE)
        {   
            RmgTimer RTL("6-Davidson: ldaUop x psi"); 
            newsint = this->orbitalsint_local + nbase * this->OrbitalProjector->get_num_nonloc_ions() * 
                this->OrbitalProjector->get_pstride() * ct.noncoll_factor;
            LdaplusUxpsi(this, nbase, notconv, newsint);
        }
        RT1 = new RmgTimer("6-Davidson: apply hamiltonian");
        ApplyHamiltonianBlock<KpointType> (this, nbase, notconv, h_psi, vtot, vxc_psi);
        delete RT1;


        // Update the reduced Hamiltonian and S matrices
        RT1 = new RmgTimer("6-Davidson: matrix setup/reduce");
        RmgGemm(trans_a, trans_n, nbase+notconv, notconv, pbasis_noncoll, alphavel, psi, pbasis_noncoll, &h_psi[nbase*pbasis_noncoll], pbasis_noncoll, beta, &hr[nbase*ct.max_states], ct.max_states);

#if HAVE_ASYNC_ALLREDUCE
        // Asynchronously reduce it
        MPI_Request MPI_reqAij;
        if(ct.use_async_allreduce)
            MPI_Iallreduce(MPI_IN_PLACE, (double *)&hr[nbase*ct.max_states], notconv * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqAij);
        else
            BlockAllreduce((double *)&hr[nbase*ct.max_states], (size_t)notconv*(size_t)ct.max_states * (size_t)factor, pct.grid_comm);
#else
        BlockAllreduce((double *)&hr[nbase*ct.max_states], (size_t)notconv*(size_t)ct.max_states * (size_t)factor, pct.grid_comm);
#endif

        RmgGemm(trans_a, trans_n, nbase+notconv, notconv, pbasis_noncoll, alphavel, psi, pbasis_noncoll, &s_psi[nbase*pbasis_noncoll], pbasis_noncoll, beta, &sr[nbase*ct.max_states], ct.max_states);

#if HAVE_ASYNC_ALLREDUCE
        // Wait for Aij request to finish
        if(ct.use_async_allreduce) MPI_Wait(&MPI_reqAij, MPI_STATUS_IGNORE);
#endif

#if HAVE_ASYNC_ALLREDUCE
        // Asynchronously reduce Sij request
        MPI_Request MPI_reqSij;
        if(ct.use_async_allreduce)
           MPI_Iallreduce(MPI_IN_PLACE, (double *)&sr[nbase*ct.max_states], notconv * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqSij);
        else
            BlockAllreduce((double *)&sr[nbase*ct.max_states], (size_t)notconv*(size_t)ct.max_states * (size_t)factor, pct.grid_comm);
#else
        BlockAllreduce((double *)&sr[nbase*ct.max_states], (size_t)notconv*(size_t)ct.max_states * (size_t)factor, pct.grid_comm);
#endif

#if HAVE_ASYNC_ALLREDUCE
        // Wait for S request to finish
        if(ct.use_async_allreduce) MPI_Wait(&MPI_reqSij, MPI_STATUS_IGNORE);
#endif
        delete RT1;

        nbase = nbase + notconv;
        std::complex<double> *hr_C, *sr_C;
        hr_C = (std::complex<double> *)hr;
        sr_C = (std::complex<double> *)sr;

        for(int i=0;i < nbase;i++) {
            for(int j=i+1;j < nbase;j++) {

                if(typeid(KpointType) == typeid(std::complex<double>))
                {
                    hr_C[j + i*ct.max_states] = std::conj(hr_C[i + j*ct.max_states]);
                    sr_C[j + i*ct.max_states] = std::conj(sr_C[i + j*ct.max_states]);
                }
                else
                {
                    hr[j + i*ct.max_states] = hr[i + j*ct.max_states];
                    sr[j + i*ct.max_states] = sr[i + j*ct.max_states];
                }
            }
        }

        RT1 = new RmgTimer("6-Davidson: diagonalization");
        int info = GeneralDiag(hr, sr, eigsw, vr, nbase, nstates, ct.max_states, ct.subdiag_driver);
        delete RT1;
        if(info) {
            if(pct.gridpe == 0) printf("\n WARNING: Davidson GeneralDiag info = %d", info);
            return;
            throw RmgFatalException() << info<<" " <<nstates <<" " << nbase << " Diagonalization failed in Davidson, terminating." << " in " << __FILE__ << " at line " << __LINE__ << "\n";
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
            double occ = this->Kstates[st].occupation[0];
            if(ct.spin_flag) occ+= this->Kstates[st].occupation[1];
            tolerances[st] = fabs(eigs[st] - eigsw[st]);
            if((tolerances[st] > max_tol) && (fabs(occ) > 0.0002)) {
                max_tol = tolerances[st];
                max_tol_state = st;
                avg_occ_tol += tolerances[st];
                occ_states++;
            }
            if((tolerances[st] < min_tol) && (fabs(occ) <= 0.0002)) {
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
            double occ = this->Kstates[st].occupation[0];
            if(ct.spin_flag) occ+= this->Kstates[st].occupation[1];
            if(fabs(occ) > 0.0002) {
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
        if(((steps == (ct.david_max_steps-1)) || ((nbase+notconv) > ct.max_states) || (notconv == 0))) {

            // Rotate orbitals
            RT1 = new RmgTimer("6-Davidson: rotate orbitals");
#if CUDA_ENABLED || HIP_ENABLED
            KpointType *npsi = (KpointType *)RmgMallocHost(nstates*pbasis_noncoll*sizeof(KpointType));
#else
            KpointType *npsi = new KpointType[nstates*pbasis_noncoll];
#endif
            RmgGemm(trans_n, trans_n, pbasis_noncoll, nstates, nbase, alpha, psi, pbasis_noncoll, vr, ct.max_states, beta, npsi, pbasis_noncoll);
            for(int idx=0;idx < nstates*pbasis_noncoll;idx++)psi[idx] = npsi[idx];
#if CUDA_ENABLED || HIP_ENABLED
            RmgFreeHost(npsi);
#else
            delete [] npsi;
#endif
            delete RT1;


            if(notconv == 0) {
                // We use a single non update davidson cycle to get a variational value for the total
                // energy when the multigrid solver is used so we don't want to write any davidson info
                // in this case
                if (!Verify ("kohn_sham_solver","multigrid", this->ControlMap))
                    rmg_printf("Davidson converged in %d steps\n", steps+1);
                break;  // done
            }

            if(steps == (ct.david_max_steps-1)) {
                if (!Verify ("kohn_sham_solver","multigrid", this->ControlMap))
                    rmg_printf("Davidson incomplete convergence steps = %d\n", steps + 1);
                break;
            }

            // refresh s_psi and h_psi
            RT1 = new RmgTimer("6-Davidson: refresh h_psi and s_psi");
            RmgGemm(trans_n, trans_n, pbasis_noncoll, nstates, nbase, alpha, s_psi, pbasis_noncoll, vr, ct.max_states, beta, &psi[nstates*pbasis_noncoll], pbasis_noncoll);
            if(!ct.norm_conserving_pp) for(int idx=0;idx < nstates*pbasis_noncoll;idx++)s_psi[idx] = psi[nstates*pbasis_noncoll + idx];

            RmgGemm(trans_n, trans_n, pbasis_noncoll, nstates, nbase, alpha, h_psi, pbasis_noncoll, vr, ct.max_states, beta, &psi[nstates*pbasis_noncoll], pbasis_noncoll);
            for(int idx=0;idx < nstates*pbasis_noncoll;idx++)h_psi[idx] = psi[nstates*pbasis_noncoll + idx];
            delete RT1;


            // Reset hr,sr,vr
            RT1 = new RmgTimer("6-Davidson: reset hr,sr,vr");
            nbase = nstates;
            for(int ix=0;ix < ct.max_states*ct.max_states;ix++) hr[ix] = KpointType(0.0);
            for(int ix=0;ix < ct.max_states*ct.max_states;ix++) sr[ix] = KpointType(0.0);
            for(int ix=0;ix < ct.max_states*ct.max_states;ix++) vr[ix] = KpointType(0.0);
            for(int st=0;st < nbase;st++) {
                hr[st + st*ct.max_states] = eigs[st];
                sr[st + st*ct.max_states] = KpointType(1.0);
                vr[st + st*ct.max_states] = KpointType(1.0);
            }
            delete RT1;

        }

    }

    // Copy eigs from compact array back into state structure
    for(int st = 0;st < nstates;st++) this->Kstates[st].eig[0] = eigs[st];



#if CUDA_ENABLED || HIP_ENABLED
    RmgFreeHost(vr);
    RmgFreeHost(sr);
    RmgFreeHost(hr);
    RmgFreeHost(h_psi);
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
    this->BetaProjector->project(this, this->newsint_local, 0, nstates*ct.noncoll_factor, weight);
    delete RT1;

}

