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
#include "make_conf.h"
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "blas.h"

#include "prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#include "Elpa.h"


template char * Subdiag_Elpa<double> (Kpoint<double> *kptr, double *Aij, double *Bij, double *Sij, double *eigs, double *eigvectors);
template char * Subdiag_Elpa<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *Aij, std::complex<double> *Bij, std::complex<double> *Sij, double *eigs, std::complex<double> *eigvectors);

Elpa *MainElpa;

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
template <typename KpointType>
char * Subdiag_Elpa (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors)
{

#if !SCALAPACK_LIBS || !USE_ELPA
    rmg_printf("This version of RMG was not built with Scalapack/Elpa support. Redirecting to LAPACK.");
    return Subdiag_Lapack(kptr, Aij, Bij, Sij, eigs, eigvectors);
#else

    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);

    static char *trans_n = "n";
    static char *trans_t = "t";

//  folded spectrum with scalapack is experimental. Uncomment the second line if you want to try it.
    bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));
    //use_folded = false;

    int ione=1;
    int num_states = kptr->nstates;
    int factor = 1;
    if(!ct.is_gamma) factor=2;

    // Create 1 scalapack instance per grid_comm. We use a static Scalapack here since initialization on large systems is expensive
    if(!MainElpa) {
        // Need some code here to decide how to set the number of scalapack groups but for now use just 1
        int scalapack_groups = 1;
        int last = !ct.use_folded_spectrum;
        MainElpa = new Elpa(scalapack_groups, pct.thisimg, ct.images_per_node, num_states,
                     ct.scalapack_block_factor, last, pct.grid_comm);
        MainElpa->GetCommunicators();
    }

    bool participates = MainElpa->Participates();
    int scalapack_nprow = MainElpa->GetRows();
    int scalapack_npcol = MainElpa->GetCols();
    int scalapack_npes = scalapack_nprow * scalapack_npcol;
    int root_npes = MainElpa->GetRootNpes();
    int elpa_comm_rows = MainElpa->GetElpaCommRows();
    int elpa_comm_cols = MainElpa->GetElpaCommCols();

    // Allocate and clear distributed matrices */
    static KpointType *distAij;
    static KpointType *distBij;
    static KpointType *distSij;
    static KpointType *distCij;
    KpointType *Cij = new KpointType[num_states * num_states]();


    RmgTimer *RT1;

    int *desca;
    int dist_length;
    
    if (participates) {

        dist_length = MainElpa->GetDistMdim() * MainElpa->GetDistNdim();
        desca = MainElpa->GetDistDesca();


        if(!distAij) {
            int retval1 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distAij);
            int retval2 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distBij);
            int retval3 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distSij);
            int retval4 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distCij);
            if((retval1 != MPI_SUCCESS) || (retval2 != MPI_SUCCESS) || (retval3 != MPI_SUCCESS) || (retval4 != MPI_SUCCESS)) {
                rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in Subdiag_Scalapack");
            }
        }

        // Copy matrices to dist arrays
        RT1 = new RmgTimer("4-Diagonalization: distribute matrices.");
        MainElpa->CopySquareMatrixToDistArray(Aij, distAij, num_states, desca);
        MainElpa->CopySquareMatrixToDistArray(Sij, distSij, num_states, desca);
        MainElpa->CopySquareMatrixToDistArray(eigvectors, distBij, num_states, desca);
        delete(RT1);

        // Create unitary matrix
        for (int idx = 0; idx < num_states; idx++) {
            Cij[idx * num_states + idx] = ONE_t;
        }

        // distribute unitary matrix
        MainElpa->CopySquareMatrixToDistArray(Cij, distCij, num_states, desca);

        if(!ct.norm_conserving_pp || (ct.norm_conserving_pp && ct.discretization == MEHRSTELLEN_DISCRETIZATION)) {

            RT1 = new RmgTimer("4-Diagonalization: Invert Bij");
            // Get matrix that is inverse to B
            {
                int info=0;

                int ipiv_size = MainElpa->GetIpivSize();
                int *ipiv = new int[ipiv_size]();

                /*Inverse of B should be in Cij */
                MainElpa->Pgesv (&num_states, &num_states, distBij, &ione, &ione, desca, ipiv, distCij, &ione,
                            &ione, desca, &info);

                if (info)
                {
                    rmg_printf ("\n PE %d: p{d,z}gesv failed, info is %d", pct.gridpe, info);
                    rmg_error_handler (__FILE__, __LINE__, " p{d,z}gesv failed");
                }

                delete [] ipiv;
            }
            delete(RT1);


            RT1 = new RmgTimer("4-Diagonalization: matrix setup");
            /*Multiply inverse of B and and A */
            {
                char *trans = "n";
                KpointType alpha(1.0);
                KpointType beta(0.0);

                /*B^-1*A */
                MainElpa->Pgemm(trans, trans, &num_states, &num_states, &num_states, &alpha,
                            distCij, &ione, &ione, desca, distAij, &ione,
                            &ione, desca, &beta, distBij, &ione, &ione, desca);

                /*Multiply the result with Sij, result is in distCij */
                MainElpa->Pgemm (trans, trans, &num_states, &num_states, &num_states, &alpha,
                            distSij, &ione, &ione, desca, distBij, &ione, 
                            &ione, desca, &beta, distCij, &ione, &ione, desca);

                // Copy result into Bij and Aij
                for(int idx=0;idx < dist_length;idx++) distBij[idx] = distCij[idx];

            }
            delete(RT1);

        }
        else {

            // For norm conserving S=B so no need to invert and S*(B-1)*A=A so just copy A into Cij
            // to pass to eigensolver
            for(int ix=0;ix < dist_length;ix++) distBij[ix] = distAij[ix];

        }

    }


    if(use_folded) {
#if 0
        // We have to gather distBij back to Bij and then broadcast it to all nodes in the root
        // Sij is still present on all nodes in original form
        //MainElpa->GatherMatrix(Bij, distBij);
        //MainElpa->BcastRoot(Bij, factor * num_states * num_states, MPI_DOUBLE);

        FoldedSpectrumScalapack<double> ((Kpoint<double> *)kptr, num_states, (double *)Bij, (double *)distBij, num_states, (double *)Sij, num_states, eigs, (double *)Aij, MainElpa, SUBDIAG_LAPACK, ct.scalapack_block_factor);

        for(int idx=0;idx< num_states * num_states;idx++)eigvectors[idx] = Bij[idx];
        // Broadcast results if required
        if(root_npes != scalapack_npes) { 
            RT1 = new RmgTimer("4-Diagonalization: MPI_Bcast");
            MainElpa->BcastRoot(eigvectors, factor * num_states * num_states, MPI_DOUBLE);
            MainElpa->BcastRoot(eigs, num_states, MPI_DOUBLE);
            delete(RT1);
        }
#endif
    }
    else {

        if(participates) {

            /****************** Find Matrix of Eigenvectors *****************************/
            /* Using lwork=-1, PDSYGVX should return minimum required size for the work array */
            RT1 = new RmgTimer("4-Diagonalization: ELPA");
            {
                int info, ibtype=1;
                char *uplo = "u", *jobz = "v";

                if(ct.is_gamma) {

                    int useQr = false;
                    int useGPU = false;
                    int wantDebug = true;
                    double scale=1.0, rone = 1.0;
                    KpointType alpha(1.0);
                    KpointType beta(0.0);
                    int THIS_REAL_ELPA_KERNEL_API = ELPA2_REAL_KERNEL_SSE;

printf("DDD  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d\n",num_states,MainElpa->GetRow(),MainElpa->GetCol(),ct.scalapack_block_factor,MainElpa->GetDistMdim(),MainElpa->GetDistNdim(), elpa_comm_rows, elpa_comm_cols, MainElpa->GetRow(),MainElpa->GetCol());

//                    info = elpa_cholesky_real_double(num_states, (double *)distSij, 
//                           MainElpa->GetDistMdim(), ct.scalapack_block_factor, MainElpa->GetDistNdim(), 
//                           elpa_comm_rows, elpa_comm_cols, wantDebug);
                    pdpotrf_(uplo, &num_states, (double *)distSij,  &ione, &ione, desca,  &info);

//                    info = elpa_invert_trm_real_double(num_states, (double *)distSij, 
//                           MainElpa->GetDistMdim(), ct.scalapack_block_factor, MainElpa->GetDistNdim(),
//                           elpa_comm_rows, elpa_comm_cols, wantDebug);

                    // Get pdsyngst_ workspace
                    int lwork = -1;
                    double lwork_tmp;
                    pdsyngst_(&ibtype, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                            (double *)distSij, &ione, &ione, desca, &scale, &lwork_tmp, &lwork, &info);
                    lwork = 2*(int)lwork_tmp;
                    double *work2 = new double[lwork];

                    pdsyngst_(&ibtype, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                            (double *)distSij, &ione, &ione, desca, &scale, work2, &lwork, &info);


#if 0
                    info = elpa_mult_at_b_real_double('U', 'F', num_states, num_states,
                           (double *)distBij, MainElpa->GetDistMdim(), MainElpa->GetDistNdim(),
                           (double *)distSij, MainElpa->GetDistMdim(), MainElpa->GetDistNdim(),
                           ct.scalapack_block_factor, elpa_comm_rows, elpa_comm_cols,
                           (double *)distCij, MainElpa->GetDistMdim(), MainElpa->GetDistNdim());

                    info = elpa_mult_at_b_real_double('U', 'F', num_states, num_states,
                           (double *)distCij, MainElpa->GetDistMdim(), MainElpa->GetDistNdim(),
                           (double *)distSij, MainElpa->GetDistMdim(), MainElpa->GetDistNdim(),
                           ct.scalapack_block_factor, elpa_comm_rows, elpa_comm_cols,
                           (double *)distBij, MainElpa->GetDistMdim(), MainElpa->GetDistNdim());
                MainElpa->Pgemm (trans_n, trans_n, &num_states, &num_states, &num_states, &alpha,
                            distSij, &ione, &ione, desca, distBij, &ione, 
                            &ione, desca, &beta, distCij, &ione, &ione, desca);
                MainElpa->Pgemm (trans_n, trans_t, &num_states, &num_states, &num_states, &alpha,
                            distCij, &ione, &ione, desca, distSij, &ione, 
                            &ione, desca, &beta, distBij, &ione, &ione, desca);
#endif
                    // Get workspace required
                    lwork = -1;
                    int liwork=-1;
                    int liwork_tmp;
                    pdsyevd_(jobz, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                            eigs, (double *)distAij, &ione, &ione, desca, &lwork_tmp, &lwork, &liwork_tmp, &liwork, &info);
                    lwork = 16*(int)lwork_tmp;
                    liwork = 16*num_states;
                    double *nwork = new double[lwork];
                    int *iwork = new int[liwork];

                    // and now solve it 
//                    pdsyevd_(jobz, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
//                            eigs, (double *)distAij, &ione, &ione, desca, nwork, &lwork, iwork, &liwork, &info);
printf("DDDD  %d  %d  %d  %d  %d  %d  %d  %d  %d\n",desca[0],desca[1],desca[2],desca[3],desca[4],desca[5],desca[6],desca[7],desca[8]);
printf("MMMM  %d  %d  %d  %d\n",MainElpa->GetDistMdim(),MainElpa->GetDistNdim(),num_states,ct.scalapack_block_factor);
                    info = elpa_solve_evp_real_2stage_double_precision(
                           num_states, num_states, 
                           (double *)distBij, MainElpa->GetDistMdim(),
                           eigs, 
                           (double *)distAij, MainElpa->GetDistMdim(),
                           ct.scalapack_block_factor, MainElpa->GetDistNdim(), 
                           elpa_comm_rows, elpa_comm_cols, MPI_Comm_c2f(MainElpa->GetRootComm()), 
                           THIS_REAL_ELPA_KERNEL_API, useQr, useGPU);

for(int i=0;i<num_states;i++){
printf("EIGS = %f   %d\n",Ha_eV*eigs[i],info);
}
printf("INFO4 = %d\n",info);

                }
                else {

                }


            }
            delete(RT1);
            
            // Clear eigvectors
            for (int idx = 0; idx < num_states*num_states; idx++) {
                eigvectors[idx] = ZERO_t;
            }
            // Gather distributed results from distAij into eigvectors
            //MainElpa->GatherMatrix(eigvectors, distAij);
            MainElpa->CopyDistArrayToSquareMatrix(eigvectors, distAij, num_states, desca);
            MainElpa->Allreduce(MPI_IN_PLACE, eigvectors, factor *num_states*num_states, MPI_DOUBLE, MPI_SUM);


        }

        // Finally send eigenvalues and vectors to everyone 
        RT1 = new RmgTimer("4-Diagonalization: MPI_Bcast");
        MainElpa->BcastRoot(eigvectors, factor * num_states * num_states, MPI_DOUBLE);
        MainElpa->BcastRoot(eigs, num_states, MPI_DOUBLE);
        delete(RT1);

    }


    delete [] Cij;

//    if(use_folded) return trans_t;   // Currently using pdsyngst in lower level routine. If
//    switch to FOLDED_GSE must uncomment
    return trans_n;
#endif

}
