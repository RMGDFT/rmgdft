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

#include "Scalapack.h"


template char * Subdiag_Scalapack<double> (Kpoint<double> *kptr, double *Aij, double *Bij, double *Sij, double *eigs, double *eigvectors);
template char * Subdiag_Scalapack<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *Aij, std::complex<double> *Bij, std::complex<double> *Sij, double *eigs, std::complex<double> *eigvectors);

Scalapack *MainSp;

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
template <typename KpointType>
char * Subdiag_Scalapack (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors)
{

#if !SCALAPACK_LIBS
    rmg_printf("This version of RMG was not built with Scalapack support. Redirecting to LAPACK.");
    return Subdiag_Lapack(kptr, Aij, Bij, Sij, eigs, eigvectors);
#else

    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);

    //static char *trans_t = "t";
    static char *trans_n = "n";

//  folded spectrum with scalapack is experimental. Uncomment the second line if you want to try it.
    bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));
    //use_folded = false;

    int ione=1;
    int num_states = kptr->nstates;
    int factor = 1;
    if(!ct.is_gamma) factor=2;

    // Create 1 scalapack instance per grid_comm. We use a static Scalapack here since initialization on large systems is expensive
    if(!MainSp) {
        // Need some code here to decide how to set the number of scalapack groups but for now use just 1
        int scalapack_groups = 1;
        int last = !ct.use_folded_spectrum;
        MainSp = new Scalapack(scalapack_groups, pct.thisimg, ct.images_per_node, num_states,
                     ct.scalapack_block_factor, last, pct.grid_comm);

    }

    bool participates = MainSp->Participates();
    int scalapack_nprow = MainSp->GetRows();
    int scalapack_npcol = MainSp->GetCols();
    int scalapack_npes = scalapack_nprow * scalapack_npcol;
    int root_npes = MainSp->GetRootNpes();

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

        dist_length = MainSp->GetDistMdim() * MainSp->GetDistNdim();
        desca = MainSp->GetDistDesca();


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
        MainSp->CopySquareMatrixToDistArray(Aij, distAij, num_states, desca);
        MainSp->CopySquareMatrixToDistArray(Sij, distSij, num_states, desca);
        MainSp->CopySquareMatrixToDistArray(eigvectors, distBij, num_states, desca);
        delete(RT1);

        // Create unitary matrix
        for (int idx = 0; idx < num_states; idx++) {
            Cij[idx * num_states + idx] = ONE_t;
        }

        // distribute unitary matrix
        MainSp->CopySquareMatrixToDistArray(Cij, distCij, num_states, desca);

        if(!ct.norm_conserving_pp || (ct.norm_conserving_pp && ct.discretization == MEHRSTELLEN_DISCRETIZATION)) {

            RT1 = new RmgTimer("4-Diagonalization: Invert Bij");
            // Get matrix that is inverse to B
            {
                int info=0;

                int ipiv_size = MainSp->GetIpivSize();
                int *ipiv = new int[ipiv_size]();

                /*Inverse of B should be in Cij */
                MainSp->Pgesv (&num_states, &num_states, distBij, &ione, &ione, desca, ipiv, distCij, &ione,
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
                MainSp->Pgemm(trans, trans, &num_states, &num_states, &num_states, &alpha,
                            distCij, &ione, &ione, desca, distAij, &ione, 
                            &ione, desca, &beta, distBij, &ione, &ione, desca);

                /*Multiply the result with Sij, result is in distCij */
                MainSp->Pgemm (trans, trans, &num_states, &num_states, &num_states, &alpha,
                            distSij, &ione, &ione, desca, distBij, &ione, 
                            &ione, desca, &beta, distCij, &ione, &ione, desca);

                // Copy result into Bij
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

        // We have to gather distBij back to Bij and then broadcast it to all nodes in the root
        // Sij is still present on all nodes in original form
        //MainSp->GatherMatrix(Bij, distBij);
        //MainSp->BcastRoot(Bij, factor * num_states * num_states, MPI_DOUBLE);

        FoldedSpectrumScalapack<double> ((Kpoint<double> *)kptr, num_states, (double *)Bij, (double *)distBij, num_states, (double *)Sij, num_states, eigs, (double *)Aij, MainSp, SUBDIAG_LAPACK, ct.scalapack_block_factor);

        for(int idx=0;idx< num_states * num_states;idx++)eigvectors[idx] = Bij[idx];
        // Broadcast results if required
        if(root_npes != scalapack_npes) { 
            RT1 = new RmgTimer("4-Diagonalization: MPI_Bcast");
            MainSp->BcastRoot(eigvectors, factor * num_states * num_states, MPI_DOUBLE);
            MainSp->BcastRoot(eigs, num_states, MPI_DOUBLE);
            delete(RT1);
        }

    }
    else {

        if(participates) {

            /****************** Find Matrix of Eigenvectors *****************************/
            /* Using lwork=-1, PDSYGVX should return minimum required size for the work array */
            RT1 = new RmgTimer("4-Diagonalization: PDSYGVX/PZHEGVX");
            {
                char *uplo = "l", *jobz = "v";
                int info;

                if(ct.is_gamma) {

                    int ibtype = 1;
                    double scale=1.0, rone = 1.0;

                    pdpotrf_(uplo, &num_states, (double *)distSij,  &ione, &ione, desca,  &info);

                    // Get pdsyngst_ workspace
                    int lwork = -1;
                    double lwork_tmp;
                    pdsyngst_(&ibtype, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                            (double *)distSij, &ione, &ione, desca, &scale, &lwork_tmp, &lwork, &info);
                    lwork = 2*(int)lwork_tmp; 
                    double *work2 = new double[lwork];
                    
                    pdsyngst_(&ibtype, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                            (double *)distSij, &ione, &ione, desca, &scale, work2, &lwork, &info);

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
                    pdsyevd_(jobz, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                            eigs, (double *)distAij, &ione, &ione, desca, nwork, &lwork, iwork, &liwork, &info);

                    pdtrsm_("Left", uplo, "T", "N", &num_states, &num_states, &rone, (double *)distSij, &ione, &ione, desca,
                            (double *)distAij, &ione, &ione, desca);
                    delete [] iwork;
                    delete [] nwork;
                    delete [] work2;

                }
                else {

                    int ibtype = 1;
                    double scale=1.0, rone[2] = {1.0, 0.0};

                    pzpotrf_(uplo, &num_states, (double *)distSij,  &ione, &ione, desca,  &info);

                    pzhegst_(&ibtype, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                                (double *)distSij, &ione, &ione, desca, &scale, &info);

                    // Get workspace required
                    int lwork = -1, liwork=-1, lrwork=-1;
                    double lwork_tmp[2], lrwork_tmp;
                    int liwork_tmp;
                    pzheevd_(jobz, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                                eigs, (double *)distAij, &ione, &ione, desca, lwork_tmp, &lwork, &lrwork_tmp, &lrwork, &liwork_tmp, &liwork, &info);
                    lwork = (int)lwork_tmp[0]+1;
                    liwork = 16*num_states;
                    lrwork = 2*(int)lrwork_tmp;
                    double *rwork = new double[lrwork];
                    double *nwork = new double[lwork*2];
                    int *iwork = new int[liwork];

                    // and now solve it
                    pzheevd_(jobz, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                                eigs, (double *)distAij, &ione, &ione, desca, nwork, &lwork, (double *)rwork, &lrwork, iwork, &liwork, &info);

                    pztrsm_("Left", uplo, "C", "N", &num_states, &num_states, rone, (double *)distSij, &ione, &ione, desca,
                                (double *)distAij, &ione, &ione, desca);

                    delete [] iwork;
                    delete [] nwork;
                    delete [] rwork;

                }

                if (info)
                {
                    rmg_printf ("\n PDSYGVX failed, info is %d", info);
                    rmg_error_handler (__FILE__, __LINE__, "PDSYGVX failed");
                }


            }
            delete(RT1);
            
            // Clear eigvectors
            for (int idx = 0; idx < num_states*num_states; idx++) {
                eigvectors[idx] = ZERO_t;
            }
            // Gather distributed results from distAij into eigvectors
            //MainSp->GatherMatrix(eigvectors, distAij);
            MainSp->CopyDistArrayToSquareMatrix(eigvectors, distAij, num_states, desca);
            MainSp->Allreduce(MPI_IN_PLACE, eigvectors, factor *num_states*num_states, MPI_DOUBLE, MPI_SUM);


        }

        // Finally send eigenvalues and vectors to everyone 
        RT1 = new RmgTimer("4-Diagonalization: MPI_Bcast");
        MainSp->BcastRoot(eigvectors, factor * num_states * num_states, MPI_DOUBLE);
        MainSp->BcastRoot(eigs, num_states, MPI_DOUBLE);
        delete(RT1);

    }


    delete [] Cij;

//    if(use_folded) return trans_t;   // Currently using pdsyngst in lower level routine. If
//    switch to FOLDED_GSE must uncomment
    return trans_n;
#endif

}
