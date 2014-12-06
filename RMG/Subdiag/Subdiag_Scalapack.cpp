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

    static char *trans_t = "t";
    static char *trans_n = "n";

    bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));

    int izero = 0;
    int ione=1;
    int dist_length=1;
    int num_states = kptr->nstates;
    int factor = 1;
    if(!ct.is_gamma) factor=2;

    // Create 1 scalapack instance per grid_comm. We use a static Scalapack here since initialization on large systems is expensive
    static Scalapack *MainSp;
    if(!MainSp) {
        MainSp = new Scalapack(4, pct.thisimg, ct.images_per_node, num_states, num_states, ct.scalapack_block_factor, ct.scalapack_block_factor, pct.grid_comm);

    }


    dist_length = MainSp->GetDistMdim() * MainSp->GetDistNdim();
    int *desca = MainSp->GetDistDesca();
    bool participates = MainSp->Participates();
    MPI_Comm scalapack_comm = MainSp->GetComm();
    int scalapack_nprow = MainSp->GetRows();
    int scalapack_npcol = MainSp->GetCols();

    if(dist_length == 0) dist_length = 1;   // Just to keep allocations from complaining


    // Allocate and clear distributed matrices */
    KpointType *distAij = new KpointType[dist_length]();
    KpointType *distBij = new KpointType[dist_length]();
    KpointType *distSij = new KpointType[dist_length]();
    KpointType *distCij = new KpointType[dist_length]();



    // Reduce and distribute matrices
    RmgTimer *RT1 = new RmgTimer("Diagonalization: distribute matrices.");
    MainSp->DistributeMatrix(Aij, distAij, num_states, num_states);
    MainSp->DistributeMatrix(Sij, distSij, num_states, num_states);
    MainSp->DistributeMatrix(eigvectors, distBij, num_states, num_states);
    delete(RT1);

    // Create unitary matrix
    KpointType *Cij = new KpointType[num_states * num_states]();
    for (int idx = 0; idx < num_states; idx++) {
        Cij[idx * num_states + idx] = ONE_t;
    }


    if (participates) {

        // distribute unitary matrix
        MainSp->DistributeMatrix(Cij, distCij, num_states, num_states);

        if(!ct.norm_conserving_pp || (ct.norm_conserving_pp && ct.discretization == MEHRSTELLEN_DISCRETIZATION)) {

            RT1 = new RmgTimer("Diagonalization: Invert Bij");
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


            RT1 = new RmgTimer("Diagonalization: matrix setup");
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


        /****************** Find Matrix of Eigenvectors *****************************/
        /* Using lwork=-1, PDSYGVX should return minimum required size for the work array */
        RT1 = new RmgTimer("Diagonalization: PDSYGVX/PZHEGVX");
        {
            char *range = "a";
            char *uplo = "l", *jobz = "v";
            double vx = 0.0;
            double tol = 0.0;
            int eigs_found, eigvs_found;
            double orfac = 0.0;
            int *iwork;
            double lwork_tmp, rwork_tmp, *work2, *rwork2;
            int liwork_tmp;
            int lrwork=-1;
            int info=0;

            int *ifail = new int[num_states];
            int *iclustr = new int[2 * pct.scalapack_nprow * pct.scalapack_npcol];
            double *gap = new double[pct.scalapack_nprow * pct.scalapack_npcol];
            int lwork = -1;
            int liwork = -1;

            if(ct.is_gamma) {
                PDSYGVX (&ione, jobz, range, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                         (double *)distSij, &ione, &ione, desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
                         &eigvs_found, eigs, &orfac, (double *)distAij, &ione, &ione, desca, &lwork_tmp, &lwork,
                         &liwork_tmp, &liwork, ifail, iclustr, gap, &info);
            }
            else {
                PZHEGVX (&ione, jobz, range, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                         (double *)distSij, &ione, &ione, desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
                         &eigvs_found, eigs, &orfac, (double *)distAij, &ione, &ione, desca, &lwork_tmp, &lwork,
                         &rwork_tmp, &lrwork,
                         &liwork_tmp, &liwork, ifail, iclustr, gap, &info);
            }

            if (info)
            {
                rmg_printf ("\n PDSYGVX query failed, info is %d", info);
                rmg_error_handler (__FILE__, __LINE__, "PDSYGVX query failed");
            }

            /*set lwork and liwork */
            lwork = (int) lwork_tmp + 1;
            liwork = liwork_tmp;

            work2 = new double[2*lwork];
            iwork = new int[liwork];

            tol = 1e-15;

            if(ct.is_gamma) {
                PDSYGVX (&ione, jobz, range, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                         (double *)distSij, &ione, &ione, desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
                         &eigvs_found, eigs, &orfac, (double *)distAij, &ione, &ione, desca, work2, &lwork, iwork,
                         &liwork, ifail, iclustr, gap, &info);
            }
            else {
                // use Aij for workspace
                rwork2 = (double *)Aij;
                lrwork = (int)rwork_tmp + 1;
                std::complex<double> *rwork2 = new  std::complex<double>[lrwork];
                PZHEGVX (&ione, jobz, range, uplo, &num_states, (double *)distBij, &ione, &ione, desca,
                         (double *)distSij, &ione, &ione, desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
                         &eigvs_found, eigs, &orfac, (double *)distAij, &ione, &ione, desca, work2, &lwork,
                         (double *)rwork2, &lrwork,
                         iwork, &liwork, ifail, iclustr, gap, &info);
                delete [] rwork2;
            }

            if (info)
            {
                rmg_printf ("\n PDSYGVX failed, info is %d", info);
                rmg_error_handler (__FILE__, __LINE__, "PDSYGVX failed");
            }


            delete [] ifail;
            delete [] iclustr;
            delete [] gap;
            delete [] work2;
            delete [] iwork;

        }
        delete(RT1);
        
        // Clear eigvectors
        for (int idx = 0; idx < num_states*num_states; idx++) {
            eigvectors[idx] = ZERO_t;
        }
        // Gather distributed results from distAij into eigvectors
        MainSp->GatherMatrix(eigvectors, distAij, num_states, num_states);



    }
    else {
        // Non-participating nodes need eigvectors set to zero as well
        for (int idx = 0; idx < num_states*num_states; idx++) {
            eigvectors[idx] = ZERO_t;
        }
    }

    // Finally, sum eigvectors over all PE's in this scalapack instance
    RT1 = new RmgTimer("Diagonalization: MPI_Allreduce");
    MainSp->Allreduce(MPI_IN_PLACE, eigvectors, factor * num_states * num_states, MPI_DOUBLE, MPI_SUM);
    delete(RT1);



    /*If some processors did not participate in Scalapack,
     * broadcast eigenvalues, since only Scalapack processors have updated eigenvalues*/
    if ((pct.scalapack_nprow * pct.scalapack_npcol != pct.scalapack_npes) && (ct.diag == 1))
    {
        int item;
        item = pct.thisimg % ct.images_per_node;
        item = item * scalapack_nprow * scalapack_npcol;

        int ppp;
        MPI_Comm_size(scalapack_comm, &ppp);

        MPI_Bcast (eigs, num_states, MPI_DOUBLE, item, scalapack_comm);
    }



    delete [] Cij;
    delete [] distCij;
    delete [] distSij;
    delete [] distBij;
    delete [] distAij;
    return trans_n;
#endif

}
