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
#include "RmgException.h"
#include "Subdiag.h"
#include "Solvers.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "RmgParallelFft.h"

#include "transition.h"
#include "blas.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


extern double *rho;

// Davidson diagonalization solver

template void Davidson<double>(Kpoint<double> *, double *, int &);
template void Davidson<std::complex<double> >(Kpoint<std::complex<double>> *, double *, int &);

extern "C" void dsygvd_(int *itype, char *jobz, char *uplo, int *n, double *a, int *lda, double *b, int *ldb, double *eigs, double *work, int *lwork, int *iwork, int *liwork, int *info);
int local_diag(double *hr, double *sr, double *eigs, double *vr, int nbase, int nstates);

template <typename OrbitalType>
void Davidson (Kpoint<OrbitalType> *kptr, double *vtot, int &notconv)
{
    RmgTimer RT0("Davidson"), *RT1;
    OrbitalType alpha(1.0);
    OrbitalType beta(0.0);
    double d0 = 1.0;
    if(ct.rms < 0.1) d0 = fabs(log10(ct.rms));
    double occupied_tol = std::min(ct.rms/d0, 1.0e-6);
    // Need this since the eigensolver may become unstable for very small residuals
    occupied_tol = std::max(occupied_tol, 5.0e-12);
    double unoccupied_tol = std::max( ( occupied_tol * 5.0 ), 1.0e-5 );
    if(pct.gridpe == 0)printf("OCCUPIED TOLERANCE = %20.12e\n",occupied_tol);

    RT1 = new RmgTimer("Davidson: diagonals");
    Diagonals(kptr);
    delete RT1;

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


    // For MPI routines
    int factor = 2;
    if(ct.is_gamma) factor = 1;
    OrbitalType *NULLptr = NULL;

    double *eigs = new double[ct.max_states];
    double *eigsw = new double[ct.max_states];
    bool *converged = new bool[ct.max_states]();

    OrbitalType *hr = new OrbitalType[ct.max_states * ct.max_states]();
    OrbitalType *sr = new OrbitalType[ct.max_states * ct.max_states]();
    OrbitalType *vr = new OrbitalType[ct.max_states * ct.max_states]();
    OrbitalType *ke = new OrbitalType[ct.max_states * pbasis]();
    for(int idx = 0;idx < nstates;idx++) vr[idx*ct.max_states + idx] = OrbitalType(1.0);

#if GPU_ENABLED
    cublasStatus_t custat;
    OrbitalType *h_psi = (OrbitalType *)GpuMallocHost(pbasis * ct.max_states * sizeof(OrbitalType));
#else
    OrbitalType *h_psi = new OrbitalType[pbasis * ct.max_states];
#endif

    // short version
    OrbitalType *psi = kptr->orbital_storage;

    // Apply Hamiltonian to current set of eigenvectors. At the current time
    // kptr->ns is also computed in ApplyHamiltonianBlock by AppNls and stored in kptr->ns
    RT1 = new RmgTimer("Davidson: Betaxpsi");
    Betaxpsi (kptr);
    delete RT1;
    kptr->mix_betaxpsi(0);
    RT1 = new RmgTimer("Davidson: apply hamiltonian");
    double fd_diag = ApplyHamiltonianBlock (kptr, 0, nstates, h_psi, vtot, ke); 
    delete RT1;
    OrbitalType *s_psi = kptr->ns;

    // Copy current eigs into compact array
    for(int st1 = 0;st1 < nstates;st1++) eigs[st1] = kptr->Kstates[st1].eig[0];

    // Compute A matrix
    RT1 = new RmgTimer("Davidson: matrix setup/reduce");
    RmgGemm(trans_a, trans_n, nbase, nbase, pbasis, alphavel, psi, pbasis, h_psi, pbasis, beta, hr, ct.max_states, 
            NULLptr, NULLptr, NULLptr, false, true, false, true);

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce it
    MPI_Request MPI_reqAij;
    MPI_Iallreduce(MPI_IN_PLACE, (double *)hr, nbase * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqAij);
#else
    MPI_Allreduce(MPI_IN_PLACE, (double *)hr, nbase * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
#endif

    // Compute S matrix
    RmgGemm (trans_a, trans_n, nbase, nbase, pbasis, alphavel, psi, pbasis, s_psi, pbasis, beta, sr, ct.max_states, 
             NULLptr, NULLptr, NULLptr, false, true, false, true);

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

    local_diag((double *)hr, (double *)sr, eigs, (double *)vr, nstates, nstates);
    for(int st=0;st < nstates;st++)eigsw[st] = eigs[st];

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
        if(steps == 50) {

            for(int st1=0;st1 < notconv;st1++) {
                for(int idx = 0;idx < pbasis;idx++) {
                    psi[(st1 + nbase)*pbasis + idx] = h_psi[st1*pbasis + idx] - eigsw[nbase + st1]*s_psi[st1*pbasis + idx];
                }
            }

        }
        else {

            RmgGemm(trans_n, trans_n, pbasis, notconv, nbase, alpha, s_psi, pbasis, vr, ct.max_states, beta, &psi[nbase*pbasis], pbasis, 
                    NULLptr, NULLptr, NULLptr, false, true, false, true);

            for(int st1=0;st1 < notconv;st1++) {
                for(int idx=0;idx < pbasis;idx++) psi[(st1 + nbase)*pbasis + idx] = -eigsw[nbase + st1] * psi[(st1 + nbase)*pbasis + idx];
            }

            RmgGemm(trans_n, trans_n, pbasis, notconv, nbase, alpha, h_psi, pbasis, vr, ct.max_states, alpha, &psi[nbase*pbasis], pbasis, 
                    NULLptr, NULLptr, NULLptr, false, true, false, true);

        }

        // Apply preconditioner
        for(int st1=0;st1 < notconv;st1++) {

            for(int idx = 0;idx < pbasis;idx++) {
                double scale;
scale = std::real(ke[st1*pbasis + idx] + vtot[idx] + kptr->vnl_diag[idx] - eigsw[st1+nbase]*kptr->s_diag[idx]);
scale = std::real(fd_diag + vtot[idx] + kptr->vnl_diag[idx] - eigsw[st1+nbase]*kptr->s_diag[idx]);
scale = std::real(ke[st1*pbasis + idx]);
scale = std::copysign(fabs(scale) + 0.00001, scale);
                scale = -1.0 / scale;
                double d1 = std::real(psi[(st1 + nbase)*pbasis + idx]);
                psi[(st1 + nbase)*pbasis + idx] = scale*d1;
            }
        }



#if 1
        // Normalize correction vectors
        RT1 = new RmgTimer("Davidson: normalization");
//for(int st1=0;st1<notconv;st1++)kptr->Kstates[nbase+st1].normalize(kptr->Kstates[nbase+st1].psi, st1);
        double *norms = new double[notconv]();
        for(int st1=0;st1 < notconv;st1++) {
            for(int idx=0;idx < pbasis;idx++) norms[st1] += vel * std::norm(psi[(st1 + nbase)*pbasis + idx]);
        }

        MPI_Allreduce(MPI_IN_PLACE, (double *)norms, notconv, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

        for(int st1=0;st1 < notconv;st1++) {
             norms[st1] = 1.0 / sqrt(norms[st1]);
             for(int idx=0;idx < pbasis;idx++) psi[(st1 + nbase)*pbasis + idx] *= norms[st1];
        }
        delete [] norms;
        delete RT1;
#endif

        // Apply Hamiltonian to the new vectors
        RT1 = new RmgTimer("Davidson: Betaxpsi");
        BetaxpsiPartial (kptr, nbase, notconv);
        delete RT1;
        kptr->mix_betaxpsi(0);

        RT1 = new RmgTimer("Davidson: apply hamiltonian");
        ApplyHamiltonianBlock (kptr, nbase, notconv, h_psi, vtot, ke);
        delete RT1;


        // Update the reduced Hamiltonian and S matrices
        RT1 = new RmgTimer("Davidson: matrix setup/reduce");
        RmgGemm(trans_a, trans_n, nbase+notconv, notconv, pbasis, alphavel, psi, pbasis, &h_psi[nbase*pbasis], pbasis, beta, &hr[nbase*ct.max_states], ct.max_states, 
                NULLptr, NULLptr, NULLptr, false, true, false, true);

#if HAVE_ASYNC_ALLREDUCE
        // Asynchronously reduce it
        MPI_Request MPI_reqAij;
        MPI_Iallreduce(MPI_IN_PLACE, (double *)&hr[nbase*ct.max_states], notconv * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqAij);
#else
        MPI_Allreduce(MPI_IN_PLACE, (double *)&hr[nbase*ct.max_states], notconv * ct.max_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
#endif

        RmgGemm(trans_a, trans_n, nbase+notconv, notconv, pbasis, alphavel, psi, pbasis, &s_psi[nbase*pbasis], pbasis, beta, &sr[nbase*ct.max_states], ct.max_states, 
                 NULLptr, NULLptr, NULLptr, false, true, false, true);

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

        int info = 0;
        if(pct.is_local_master) {

            RT1 = new RmgTimer("Davidson: diagonalization");
            // Increase the resources available to this proc since the others on the local node
            // will be idle
            int nthreads = ct.THREADS_PER_NODE;
            if(pct.procs_per_host > 1) nthreads = pct.ncpus;
            omp_set_num_threads(nthreads);


            {
                int *ifail = new int[nbase];
                int lwork = 6 * nbase * nbase + 6 * nbase + 2;
                int liwork = 6*nbase;
                int eigs_found;
                double *work2 = new double[2*lwork];
                int *iwork = new int[liwork];
                double vx = 0.0;
                double tol = 1.0e-15;
                int itype = 1, ione = 1;

                OrbitalType *hsave = new OrbitalType[ct.max_states*ct.max_states];
                OrbitalType *ssave = new OrbitalType[ct.max_states*ct.max_states];
                for(int i=0;i<ct.max_states*ct.max_states;i++){
                    hsave[i] = hr[i];
                    ssave[i] = sr[i];
                }

                dsygvx (&itype, "V", "I", "L", &nbase, (double *)hr, &ct.max_states, (double *)sr, &ct.max_states,
                                        &vx, &vx, &ione, &nstates,  &tol, &eigs_found, eigsw, (double *)vr, &ct.max_states, work2,
                                        &lwork, iwork, ifail, &info);


                for(int i=0;i<ct.max_states*ct.max_states;i++){
                    hr[i] = hsave[i];
                    sr[i] = ssave[i];
                }
                delete [] ssave;
                delete [] hsave;

            }

            // Reset omp_num_threads
            omp_set_num_threads(ct.THREADS_PER_NODE);

            delete RT1;
        } // end if pct.is_local_master


        // If only one proc on this host participated broadcast results to the rest
        if(pct.procs_per_host > 1) {
            int factor = 2;
            if(ct.is_gamma) factor = 1;
            MPI_Bcast(vr, factor * nstates*ct.max_states, MPI_DOUBLE, 0, pct.local_comm);
            MPI_Bcast(eigsw, nstates, MPI_DOUBLE, 0, pct.local_comm);
            MPI_Bcast(&info, 1, MPI_INT, 0, pct.local_comm);
        }
        if(info) {
            throw RmgFatalException() << "No eigenvalues found in Davidson, terminating." << " in " << __FILE__ << " at line " << __LINE__ << "\n";
        }



        // Check convergence
        int tnotconv = nstates;
        for(int st=0;st < nstates;st++) {
            //printf("EIGS = %20.12f  %20.12f\n",eigs[st], eigsw[st]);
            if(kptr->Kstates[st].is_occupied()) {
                converged[st] = (fabs(eigs[st] - eigsw[st]) < occupied_tol);
            }
            else {
                converged[st] = (fabs(eigs[st] - eigsw[st]) < unoccupied_tol);
            }
            if(converged[st]) tnotconv--;
            if((pct.gridpe==0) && (info != 0)) printf("STATE %d e0=%20.12f e1=%20.12f  TOLERANCE = %20.12f\n",st,eigs[st],eigsw[st],fabs(eigs[st] - eigsw[st]));
            if((pct.gridpe==0)) printf("STATE %d e0=%20.12f e1=%20.12f  TOLERANCE = %20.12f\n",st,eigs[st],eigsw[st],fabs(eigs[st] - eigsw[st]));
        }

        notconv = tnotconv;
        // Copy updated eigenvalues back
        for(int st=0;st < nstates;st++) eigs[st] = eigsw[st];

        // Check if we converged to the desired tolerance and if so return. If we
        // have exceeded the maximum number of iterations then we need to do something else.
        // If the expanded basis is getting too large then we need to rotate the orbitals
        // and start the davidson iteration again.
        if(((steps == (max_steps-1)) || ((nbase+notconv) > ct.max_states) || (notconv == 0))) {

            // Rotate orbitals
            RT1 = new RmgTimer("Davidson: rotate orbitals");
            OrbitalType *npsi = new OrbitalType[nstates*pbasis];
            RmgGemm(trans_n, trans_n, pbasis, nstates, nbase, alpha, psi, pbasis, vr, ct.max_states, beta, npsi, pbasis, 
                NULLptr, NULLptr, NULLptr, false, true, false, true);
            for(int idx=0;idx < nstates*pbasis;idx++)psi[idx] = npsi[idx];
            delete [] npsi;
            delete RT1;

            if(notconv == 0) {
                if(pct.gridpe == 0) printf("Davidson converged in %d steps\n", steps+1);
                break;  // done
            }

            if(steps == (max_steps-1)) {
                if(pct.gridpe == 0) printf("Davidson incomplete convergence steps = %d\n", steps + 1);
                // Incomplete convergence, what should we do here?
                //throw RmgFatalException() << "Davidson failed to converge, terminating." << " in " << __FILE__ << " at line " << __LINE__ << "\n";

            }

            // refresh s_psi and h_psi
            RT1 = new RmgTimer("Davidson: refresh h_psi and s_psi");
            RmgGemm(trans_n, trans_n, pbasis, nstates, nbase, alpha, s_psi, pbasis, vr, ct.max_states, beta, &psi[nstates*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, true, false, true);
            for(int idx=0;idx < nstates*pbasis;idx++)s_psi[idx] = psi[nstates*pbasis + idx];

            RmgGemm(trans_n, trans_n, pbasis, nstates, nbase, alpha, h_psi, pbasis, vr, ct.max_states, beta, &psi[nstates*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, true, false, true);
            for(int idx=0;idx < nstates*pbasis;idx++)h_psi[idx] = psi[nstates*pbasis + idx];

            RmgGemm(trans_n, trans_n, pbasis, nstates, nbase, alpha, ke, pbasis, vr, ct.max_states, beta, &psi[nstates*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, true, false, true);
            for(int idx=0;idx < nstates*pbasis;idx++)ke[idx] = psi[nstates*pbasis + idx];

            delete RT1;


            // Reset hr,sr,vr
            RT1 = new RmgTimer("Davidson: reset hr,sr,vr");
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
    GpuFreeHost(h_psi);
#else
    delete [] h_psi;
#endif

    delete [] ke;
    delete [] vr;
    delete [] sr;
    delete [] hr;
    delete [] converged;
    delete [] eigsw;
    delete [] eigs;


    RT1 = new RmgTimer("Davidson: Betaxpsi");
    Betaxpsi (kptr);
    delete RT1;
    kptr->mix_betaxpsi(0);
    
}

int local_diag(double *hr, double *sr, double *eigs, double *vr, int nbase, int nstates)
{
        int info = 0;
        if(pct.is_local_master) {

            // Increase the resources available to this proc since the others on the local node
            // will be idle
            int nthreads = ct.THREADS_PER_NODE;
            if(pct.procs_per_host > 1) nthreads = pct.ncpus;
            omp_set_num_threads(nthreads);


            {
                int *ifail = new int[nbase];
                int lwork = 6 * nbase * nbase + 6 * nbase + 2;
                int liwork = 6*nbase;
                int eigs_found;
                double *work2 = new double[2*lwork];
                int *iwork = new int[liwork];
                double vx = 0.0;
                double tol = 1.0e-12;
                int itype = 1, ione = 1;

                double *hsave = new double[ct.max_states*ct.max_states];
                double *ssave = new double[ct.max_states*ct.max_states];
                for(int i=0;i<ct.max_states*ct.max_states;i++){
                    hsave[i] = hr[i];
                    ssave[i] = sr[i];
                }

//                dsygvx (&itype, "V", "I", "L", &nbase, hr, &ct.max_states, sr, &ct.max_states,
//                                        &vx, &vx, &ione, &nstates,  &tol, &eigs_found, eigs, vr, &ct.max_states, work2,
//                                        &lwork, iwork, ifail, &info);
dsygvd_(&itype, "V", "L", &nstates, hr, &ct.max_states, sr, &ct.max_states, eigs, work2, &lwork, iwork, &liwork, &info);
for(int ix=0;ix < nstates*ct.max_states;ix++)vr[ix] = hr[ix];


                for(int i=0;i<ct.max_states*ct.max_states;i++){
                    hr[i] = hsave[i];
                    sr[i] = ssave[i];
                }
                delete [] ssave;
                delete [] hsave;

            }

            // Reset omp_num_threads
            omp_set_num_threads(ct.THREADS_PER_NODE);

        } // end if pct.is_local_master


        // If only one proc on this host participated broadcast results to the rest
        if(pct.procs_per_host > 1) {
            int factor = 2;
            if(ct.is_gamma) factor = 1;
            MPI_Bcast(vr, factor * nstates*ct.max_states, MPI_DOUBLE, 0, pct.local_comm);
            MPI_Bcast(eigs, nstates, MPI_DOUBLE, 0, pct.local_comm);
            MPI_Bcast(&info, 1, MPI_INT, 0, pct.local_comm);
        }
        return info;
}
