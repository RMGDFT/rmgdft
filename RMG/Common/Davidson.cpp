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
#include "Subdiag.h"
#include "Solvers.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"

#include "transition.h"
#include "blas.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

int firstflag;
// Davidson diagonalization solver

template void Davidson<double>(Kpoint<double> *, double *);
template void Davidson<std::complex<double> >(Kpoint<std::complex<double>> *, double *);

template <typename OrbitalType>
void Davidson (Kpoint<OrbitalType> *kptr, double *vtot)
{
    RmgTimer RT0("Davidson"), *RT1;

    OrbitalType alpha(1.0);
    OrbitalType beta(0.0);
    double occupied_tol = std::min(ct.rms/100000.0, 1.0e-8);
occupied_tol = 1.0e-12;
    double unoccupied_tol = std::max( ( occupied_tol * 5.0 ), 5.0e-1 );

    int pbasis = kptr->pbasis;
    int nstates = kptr->nstates;
    int notconv = nstates;
    int nbase = nstates;
    int max_steps = 20;
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
#if 1
static OrbitalType *sdiag;
static OrbitalType *nldiag;

if(!firstflag) {
sdiag = new OrbitalType[pbasis];
nldiag = new OrbitalType[pbasis];
OrbitalType *s1 = new OrbitalType[pbasis];
for(int i=0;i<pbasis;i++)s1[i]=kptr->Kstates[0].psi[i];
for(int i=0;i<pbasis;i++)kptr->Kstates[0].psi[i] = OrbitalType(0.0);
kptr->nstates = 1;
ct.num_states = 1;

for(int i=0;i<pbasis;i++){
    kptr->Kstates[0].psi[i] = OrbitalType(1.0);
    Betaxpsi (kptr);
    AppNls(kptr, kptr->newsint_local, kptr->Kstates[0].psi, kptr->nv, kptr->ns, kptr->Bns, 0, 1);
    sdiag[i] = kptr->ns[i];
    nldiag[i] = kptr->nv[i];
    kptr->Kstates[0].psi[i] = OrbitalType(0.0);
//printf("STATUS %d\n",i);
}
kptr->nstates = nstates;
ct.num_states = kptr->nstates;

for(int i=0;i<pbasis;i++)kptr->Kstates[0].psi[i]=s1[i];
delete [] s1;
firstflag = 1;
}
#endif
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
    for(int idx = 0;idx < ct.max_states;idx++) vr[idx*ct.max_states + idx] = OrbitalType(1.0);

#if GPU_ENABLED
    cublasStatus_t custat;
    OrbitalType *h_psi = (OrbitalType *)GpuMallocHost(pbasis * ct.max_states * sizeof(OrbitalType));
#else
    OrbitalType *h_psi = new OrbitalType[pbasis * ct.max_states];
#endif

    // short verstion
    OrbitalType *psi = kptr->orbital_storage;

    // Apply Hamiltonian to current set of eigenvectors. At the current time
    // kptr->ns is also computed in ApplyHamiltonianBlock by AppNls and stored in kptr->ns
    RT1 = new RmgTimer("Davidson: apply hamiltonian");
    double fd_diag = ApplyHamiltonianBlock (kptr, 0, nstates, h_psi, vtot); 
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
long idum;
idum = 1314; 
rand0 (&idum);


    for(int steps = 0;steps < max_steps;steps++) {
double time1 = my_crtc ();

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
        RmgGemm(trans_n, trans_n, pbasis, notconv, nbase, alpha, s_psi, pbasis, vr, ct.max_states, beta, &psi[nbase*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, true, false, true);

        for(int st1=0;st1 < notconv;st1++) {
            for(int idx=0;idx < pbasis;idx++) psi[(st1 + nbase)*pbasis + idx] = -eigsw[nbase + st1] * psi[(st1 + nbase)*pbasis + idx];
        }

        RmgGemm(trans_n, trans_n, pbasis, notconv, nbase, alpha, h_psi, pbasis, vr, ct.max_states, alpha, &psi[nbase*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, true, false, true);



        // Apply preconditioner
        double eps = 1.0e-4;
        for(int st1=0;st1 < notconv;st1++) {
#if 1
            for(int idx = 0;idx < pbasis;idx++) {
//                double scale = 1.0 / std::real(-(fd_diag + vtot[idx] + kptr->nl_Bweight[idx] - eigsw[st1+nbase]*(1.0+kptr->nl_weight[idx]*kptr->nl_weight[idx])));
                double scale = -1.0 / std::real((fd_diag + vtot[idx] + nldiag[idx] - eigsw[st1+nbase]*sdiag[idx]));
                if(fabs(scale) < eps) {
                    bool neg = (scale < 0.0);
                    scale = eps;
                    if(neg) scale = -scale;
                }
                psi[(st1 + nbase)*pbasis + idx] *= scale;

            }
#endif
        }
//exit(0);
#if 1
        // Normalize correction vectors
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
#endif
        // Apply Hamiltonian to the new vectors
        kptr->nstates = nbase+notconv;  // Little bit of a hack until we update Betaxpsi
        ct.num_states = kptr->nstates;

        RT1 = new RmgTimer("Davidson: Betaxpsi");
        Betaxpsi (kptr);
        delete RT1;
        kptr->mix_betaxpsi(0);

        RT1 = new RmgTimer("Davidson: apply hamiltonian");
        //ApplyHamiltonianBlock (kptr, nbase, notconv, h_psi, vtot);
        ApplyHamiltonianBlock (kptr, 0, nbase+notconv, h_psi, vtot);
        delete RT1;
        kptr->nstates = nstates;  // And back out the hack
        ct.num_states = kptr->nstates;

        // Update the reduced Hamiltonian and S matrices
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

        nbase = nbase + notconv;
        for(int i=0;i < nbase;i++) {
            for(int j=i+1;j < nbase;j++) {
                hr[j + i*ct.max_states] = hr[i + j*ct.max_states];
                sr[j + i*ct.max_states] = sr[i + j*ct.max_states];
            }
        }


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
                double tol = 1.0e-14;
                int itype = 1, info=0, ione = 1;

                OrbitalType *hsave = new OrbitalType[ct.max_states*ct.max_states];
                OrbitalType *ssave = new OrbitalType[ct.max_states*ct.max_states];
                for(int i=0;i<ct.max_states*ct.max_states;i++){
                    hsave[i] = hr[i];
                    ssave[i] = sr[i];
                }

                dsygvx (&itype, "V", "I", "L", &nbase, (double *)hr, &ct.max_states, (double *)sr, &ct.max_states,
                                        &vx, &vx, &ione, &nstates,  &tol, &eigs_found, eigsw, (double *)vr, &ct.max_states, work2,
                                        &lwork, iwork, ifail, &info);
#if 0
for(int st=0;st<nbase;st++){
  //printf(" %d  %12.8e \n ", st, vr[st*ct.max_states + st]);
  if(std::real(vr[st*ct.max_states + st]) < 0.0) {
      for(int st1=0;st1<nbase;st1++)vr[st1*ct.max_states+st] *= -1.0;
  }
}
#endif
                for(int i=0;i<ct.max_states*ct.max_states;i++){
                    hr[i] = hsave[i];
                    sr[i] = ssave[i];
                }
                delete [] ssave;
                delete [] hsave;

                printf("NBASE = %d  EIGSFOUND = %d  INFO=%d\n",nbase,eigs_found,info);
            }

            // Reset omp_num_threads
            omp_set_num_threads(ct.THREADS_PER_NODE);

        } // end if pct.is_local_master

        // If only one proc on this host participated broadcast results to the rest
        if(pct.procs_per_host > 1) {
            int factor = 2;
            if(ct.is_gamma) factor = 1;
            MPI_Bcast(vr, factor * nstates*ct.max_states, MPI_DOUBLE, 0, pct.local_comm);
            MPI_Bcast(eigsw, nstates, MPI_DOUBLE, 0, pct.local_comm);
        }


        // Copy updated eigenvalues back
        int tnotconv = nstates;
        for(int st=0;st < nstates;st++) {
           printf("EIGS = %20.12f  %20.12f\n",eigs[st], eigsw[st]);
            if(kptr->Kstates[st].is_occupied()) {
                converged[st] = (fabs(eigs[st] - eigsw[st]) < occupied_tol);
            }
            else {
                converged[st] = (fabs(eigs[st] - eigsw[st]) < unoccupied_tol);
            }
            if(converged[st]) tnotconv--;
            if(pct.gridpe==0) printf("STATE %d TOLERANCE = %20.12f\n",st,fabs(eigs[st] - eigsw[st]));
        }
        notconv = tnotconv;
        for(int st=0;st < nstates;st++) eigs[st] = eigsw[st];

        // Check if we converged to the desired tolerance and if so return. If we
        // have exceeded the maximum number of iterations then we need to do something else.
        // If the expanded basis is getting too large then we need to rotate the orbitals
        // and start the davidson iteration again.
        if(steps == (max_steps-1) || ((nbase+notconv) >= ct.max_states) || (notconv == 0)) {

            // Rotate orbitals
            RmgGemm(trans_n, trans_n, pbasis, nstates, nbase, alpha, psi, pbasis, vr, ct.max_states, beta, h_psi, pbasis, 
                NULLptr, NULLptr, NULLptr, false, true, false, true);
//for(int idx=0;idx < nstates*pbasis;idx++)printf("DIFF = %20.12f  %20.12f  %20.12f  %20.12f\n",fabs(psi[idx] - h_psi[idx]), fabs(psi[idx] + h_psi[idx]), psi[idx], h_psi[idx]);
            for(int idx=0;idx < nstates*pbasis;idx++)psi[idx] = h_psi[idx];
printf("NOTCONV = %d  %d\n",notconv,steps);
            if(notconv == 0) {
                printf("CONVERGED davidson\n");
                break;  // done
            }
//exit(0);

            if(steps == (max_steps-1)) {
                // Incomplete convergence, what should we do here?
                break;
            }

            // refresh s_psi and h_psi
            RmgGemm(trans_n, trans_n, pbasis, nstates, nbase, alpha, s_psi, pbasis, vr, ct.max_states, beta, &psi[nstates*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, true, false, true);
            for(int idx=0;idx < nstates*pbasis;idx++)s_psi[idx] = psi[nstates*pbasis + idx];

            RmgGemm(trans_n, trans_n, pbasis, nstates, nbase, alpha, h_psi, pbasis, vr, ct.max_states, beta, &psi[nstates*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, true, false, true);
            for(int idx=0;idx < nstates*pbasis;idx++)h_psi[idx] = psi[nstates*pbasis + idx];


            // Reset hr,sr,vr
            nbase = nstates;
            for(int ix=0;ix < ct.max_states*ct.max_states;ix++) hr[ix] = OrbitalType(0.0);
            for(int ix=0;ix < ct.max_states*ct.max_states;ix++) sr[ix] = OrbitalType(0.0);
            for(int ix=0;ix < ct.max_states*ct.max_states;ix++) vr[ix] = OrbitalType(0.0);
            for(int st=0;st < nstates;st++) {
                hr[st + st*ct.max_states] = eigs[st];
                sr[st + st*ct.max_states] = OrbitalType(1.0);
                vr[st + st*ct.max_states] = OrbitalType(1.0);
            }
            
        }
        rmg_printf("\n DAVIDSON STEP TIME = %10.2f\n",my_crtc () - time1);

    }

    // Copy eigs from compact array back into state structure
    for(int st = 0;st < nstates;st++) kptr->Kstates[st].eig[0] = eigs[st];



#if GPU_ENABLED
    GpuFreeHost(h_psi);
#else
    delete [] h_psi;
#endif

    delete [] vr;
    delete [] sr;
    delete [] hr;
    delete [] converged;
    delete [] eigsw;
    delete [] eigs;

    Betaxpsi (kptr);
    kptr->mix_betaxpsi(0);
    
}


