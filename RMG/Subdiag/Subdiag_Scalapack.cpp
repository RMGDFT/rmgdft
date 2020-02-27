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

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"


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

    //static char *trans_t = "t";
    static char *trans_n = "n";

//  folded spectrum with scalapack is experimental. Comment out the second line if you want to try it.
    bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));
    use_folded = false;

    RmgTimer *DiagTimer;
    static int call_count, folded_call_count;
    if(use_folded)
    {
        DiagTimer = new RmgTimer("4-Diagonalization: scalapack folded");
        folded_call_count++;
        rmg_printf("\nDiagonalization using folded scalapack for step=%d  count=%d\n\n",ct.scf_steps, folded_call_count);
    }
    else
    {
        DiagTimer = new RmgTimer("4-Diagonalization: scalapack");
        call_count++;
        rmg_printf("\nDiagonalization using scalapack for step=%d  count=%d\n\n",ct.scf_steps, call_count);
    }


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

    RmgTimer *RT1;

    int *desca;
    int dist_length;
    static int saved_dist_length;
    
    if (participates) {

        dist_length = MainSp->GetDistMdim() * MainSp->GetDistNdim();
        desca = MainSp->GetDistDesca();

        if(dist_length != saved_dist_length)
        {
            if(distSij) {MPI_Free_mem(distSij);distSij = NULL;}
            if(distBij) {MPI_Free_mem(distBij);distBij = NULL;}
            if(distAij) {MPI_Free_mem(distAij);distAij = NULL;}
        }
        if(!distAij) {
            int retval1 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distAij);
            int retval2 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distBij);
            int retval3 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distSij);
            if((retval1 != MPI_SUCCESS) || (retval2 != MPI_SUCCESS) || (retval3 != MPI_SUCCESS)) {
                rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in Subdiag_Scalapack");
            }
            saved_dist_length = dist_length;
        }

        // Copy matrices to dist arrays
        RT1 = new RmgTimer("4-Diagonalization: distribute matrices.");
        MainSp->CopySquareMatrixToDistArray(Aij, distAij, num_states, desca);
        MainSp->CopySquareMatrixToDistArray(Sij, distSij, num_states, desca);
        MainSp->CopySquareMatrixToDistArray(Bij, distBij, num_states, desca);
        delete(RT1);

        // Copy Aij into Bij to pass to eigensolver
        for(int ix=0;ix < dist_length;ix++) distBij[ix] = distAij[ix];

    }


    if(use_folded) {

        // We have to gather distBij back to Bij and then broadcast it to all nodes in the root
        // Sij is still present on all nodes in original form
        //MainSp->GatherMatrix(Bij, distBij);
        //MainSp->BcastRoot(Bij, factor * num_states * num_states, MPI_DOUBLE);

        FoldedSpectrumScalapack<double> ((Kpoint<double> *)kptr, num_states, (double *)Bij, (double *)distBij, num_states, (double *)Sij, num_states, eigs, (double *)Aij, MainSp, SUBDIAG_LAPACK, ct.scalapack_block_factor);

        size_t pstop = (size_t)num_states * (size_t)num_states;
        memcpy(eigvectors, Bij, pstop);
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

                    pdpotrf(uplo, &num_states, (double *)distSij,  &ione, &ione, desca,  &info);

                    // Get pdsyngst_ workspace
                    int lwork = -1;
                    double lwork_tmp;
                    pdsyngst(&ibtype, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                            (double *)distSij, &ione, &ione, desca, &scale, &lwork_tmp, &lwork, &info);
                    lwork = 2*(int)lwork_tmp; 
                    double *work2 = new double[lwork];
                    
                    pdsyngst(&ibtype, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                            (double *)distSij, &ione, &ione, desca, &scale, work2, &lwork, &info);

                    // Get workspace required
                    lwork = -1;
                    int liwork=-1;
                    int liwork_tmp;
                    pdsyevd(jobz, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                            eigs, (double *)distAij, &ione, &ione, desca, &lwork_tmp, &lwork, &liwork_tmp, &liwork, &info);
                    lwork = 16*(int)lwork_tmp;
                    liwork = 16*num_states;
                    double *nwork = new double[lwork];
                    int *iwork = new int[liwork];

                    // and now solve it 
                    pdsyevd(jobz, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                            eigs, (double *)distAij, &ione, &ione, desca, nwork, &lwork, iwork, &liwork, &info);

                    pdtrsm("Left", uplo, "T", "N", &num_states, &num_states, &rone, (double *)distSij, &ione, &ione, desca,
                            (double *)distAij, &ione, &ione, desca);
                    delete [] iwork;
                    delete [] nwork;
                    delete [] work2;

                }
                else {

                    int ibtype = 1;
                    double scale=1.0, rone[2] = {1.0, 0.0};

                    pzpotrf(uplo, &num_states, (double *)distSij,  &ione, &ione, desca,  &info);

                    pzhegst(&ibtype, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                                (double *)distSij, &ione, &ione, desca, &scale, &info);

                    // Get workspace required
                    int lwork = -1, liwork=-1, lrwork=-1;
                    double lwork_tmp[2], lrwork_tmp;
                    int liwork_tmp;
                    pzheevd(jobz, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                                eigs, (double *)distAij, &ione, &ione, desca, lwork_tmp, &lwork, &lrwork_tmp, &lrwork, &liwork_tmp, &liwork, &info);
                    lwork = (int)lwork_tmp[0]+1;
                    liwork = 16*num_states;
                    lrwork = 2*(int)lrwork_tmp;
                    double *rwork = new double[lrwork];
                    double *nwork = new double[lwork*2];
                    int *iwork = new int[liwork];

                    // and now solve it
                    pzheevd(jobz, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                                eigs, (double *)distAij, &ione, &ione, desca, nwork, &lwork, (double *)rwork, &lrwork, iwork, &liwork, &info);

                    pztrsm("Left", uplo, "C", "N", &num_states, &num_states, rone, (double *)distSij, &ione, &ione, desca,
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
            size_t pstop = (size_t)num_states * (size_t)num_states;
            for (size_t idx = 0; idx < pstop; idx++) {
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

    delete DiagTimer;

//    if(use_folded) return trans_t;   // Currently using pdsyngst in lower level routine. If
//    switch to FOLDED_GSE must uncomment
    return trans_n;
#endif

}
