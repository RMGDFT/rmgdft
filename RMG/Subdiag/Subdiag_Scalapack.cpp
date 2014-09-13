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



template void Subdiag_Scalapack<double> (Kpoint<double> *kptr, double *Aij, double *Bij, double *Sij, double *eigs, double *eigvectors);
template void Subdiag_Scalapack<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *Aij, std::complex<double> *Bij, std::complex<double> *Sij, double *eigs, std::complex<double> *eigvectors);

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
template <typename KpointType>
void Subdiag_Scalapack (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors)
{

#if !SCALAPACK_LIBS
    rmg_printf("This version of RMG was not built with Scalapack support. Redirecting to LAPACK.");
    Subdiag_Lapack(kptr, Aij, Bij, Sij, eigs, eigvectors);
    return;
#else

    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);
    int izero = 0;
    int ione=1;
    int dist_length=1;
    int num_states = kptr->nstates;
    int factor = 1;
    if(!ct.is_gamma) factor=2;


    if (pct.scalapack_pe)
    {

        /*Length of distributed matrices (different on each processor) */
        dist_length = NUMROC (&num_states, &pct.desca[4], &pct.scalapack_myrow, &izero,
                      &pct.scalapack_nprow) * NUMROC (&num_states, &pct.desca[4], &pct.scalapack_mycol,
                      &izero, &pct.scalapack_npcol);

        /* Every processor for which pct.scalapack_pe should have some data and so dist_length cannot be 0*/
        if (dist_length == 0)
        rmg_error_handler(__FILE__, __LINE__, " function NUMROC returned 0, that should not happen");

    }

    // Allocate and clear distributed matrices */
    KpointType *distAij = new KpointType[dist_length]();
    KpointType *distBij = new KpointType[dist_length]();
    KpointType *distSij = new KpointType[dist_length]();
    KpointType *distCij = new KpointType[dist_length]();



    // Reduce and distribute matrices
    RmgTimer *RT1 = new RmgTimer("Diagonalization: distribute matrices.");
    distribute_mat (pct.desca, (double *)Aij, (double *)distAij, &num_states);
    distribute_mat (pct.desca, (double *)Sij, (double *)distSij, &num_states);
    distribute_mat (pct.desca, (double *)eigvectors, (double *)distBij, &num_states);
    delete(RT1);

    // Create unitary matrix
    KpointType *Cij = new KpointType[num_states * num_states]();
    for (int idx = 0; idx < num_states; idx++) {
        Cij[idx * num_states + idx] = ONE_t;
    }


    if (pct.scalapack_pe) {

        // distribute unitary matrix
        distribute_mat (pct.desca, (double *)Cij, (double *)distCij, &num_states);

        if(!ct.norm_conserving_pp) {

            RT1 = new RmgTimer("Diagonalization: Invert Bij");
            // Get matrix that is inverse to B
            {
                int info=0;
                int ipiv_size = NUMROC (&pct.desca[2], &pct.desca[4], &pct.scalapack_myrow, &pct.desca[6],
                                &pct.scalapack_nprow) + pct.desca[4];
                int *ipiv = new int[ipiv_size]();

                /*Inverse of B should be in Cij */
                if(ct.is_gamma) {
                    PDGESV (&num_states, &num_states, (double *)distBij, &ione, &ione, pct.desca, ipiv, (double *)distCij, &ione,
                            &ione, pct.desca, &info);
                }
                else {
                    PZGESV (&num_states, &num_states, (double *)distBij, &ione, &ione, pct.desca, ipiv, (double *)distCij, &ione,
                            &ione, pct.desca, &info);
                }

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
                if(ct.is_gamma) {

                    PDGEMM (trans, trans, &num_states, &num_states, &num_states, (double *)&alpha,
                            (double *)distCij, &ione, &ione, pct.desca, (double *)distAij, &ione, 
                            &ione, pct.desca, (double *)&beta, (double *)distBij, &ione, &ione, pct.desca);

                    /*Multiply the result with Sij, result is in distCij */
                    PDGEMM (trans, trans, &num_states, &num_states, &num_states, (double *)&alpha,
                            (double *)distSij, &ione, &ione, pct.desca, (double *)distBij, &ione, 
                            &ione, pct.desca, (double *)&beta, (double *)distCij, &ione, &ione, pct.desca);

                }
                else {

                    PZGEMM (trans, trans, &num_states, &num_states, &num_states, (double *)&alpha,
                            (double *)distCij, &ione, &ione, pct.desca, (double *)distAij, &ione, 
                            &ione, pct.desca, (double *)&beta, (double *)distBij, &ione, &ione, pct.desca);

                    /*Multiply the result with Sij, result is in distCij */
                    PZGEMM (trans, trans, &num_states, &num_states, &num_states, (double *)&alpha,
                            (double *)distSij, &ione, &ione, pct.desca, (double *)distBij, &ione, 
                            &ione, pct.desca, (double *)&beta, (double *)distCij, &ione, &ione, pct.desca);

                }

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
        RT1 = new RmgTimer("Diagonalization: lapack");
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
                PDSYGVX (&ione, jobz, range, uplo, &num_states, (double *)distBij, &ione, &ione, pct.desca,
                         (double *)distSij, &ione, &ione, pct.desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
                         &eigvs_found, eigs, &orfac, (double *)distAij, &ione, &ione, pct.desca, &lwork_tmp, &lwork,
                         &liwork_tmp, &liwork, ifail, iclustr, gap, &info);
            }
            else {
                PZHEGVX (&ione, jobz, range, uplo, &num_states, (double *)distBij, &ione, &ione, pct.desca,
                         (double *)distSij, &ione, &ione, pct.desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
                         &eigvs_found, eigs, &orfac, (double *)distAij, &ione, &ione, pct.desca, &lwork_tmp, &lwork,
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
                PDSYGVX (&ione, jobz, range, uplo, &num_states, (double *)distBij, &ione, &ione, pct.desca,
                         (double *)distSij, &ione, &ione, pct.desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
                         &eigvs_found, eigs, &orfac, (double *)distAij, &ione, &ione, pct.desca, work2, &lwork, iwork,
                         &liwork, ifail, iclustr, gap, &info);
            }
            else {
                // use Aij for workspace
                rwork2 = (double *)Aij;
                lrwork = (int)rwork_tmp + 1;
                std::complex<double> *rwork2 = new  std::complex<double>[lrwork];
                PZHEGVX (&ione, jobz, range, uplo, &num_states, (double *)distBij, &ione, &ione, pct.desca,
                         (double *)distSij, &ione, &ione, pct.desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
                         &eigvs_found, eigs, &orfac, (double *)distAij, &ione, &ione, pct.desca, work2, &lwork,
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
        
        // Gather result onto global_matrix
        matgather ((double *)distAij, pct.desca, (double *)eigvectors, num_states);


    }
    else {
        // Non-scalapack nodes should have Aij set to zero
        for (int idx = 0; idx < num_states*num_states; idx++) {
            eigvectors[idx] = ZERO_t;
        }
    }

    // Finally, sum eigvectors over all PEs
    RT1 = new RmgTimer("Diagonalization: MPI_Allreduce");
    MPI_Allreduce(MPI_IN_PLACE, eigvectors, factor * num_states * num_states, MPI_DOUBLE, MPI_SUM, pct.scalapack_comm);
    delete(RT1);



    /*If some processors did not participate in Scalapack,
     * broadcast eigenvalues, since only Scalapack processors have updated eigenvalues*/
    if ((pct.scalapack_nprow * pct.scalapack_npcol != pct.scalapack_npes) && (ct.diag == 1))
    {
        int item;
        item = pct.thisimg % ct.images_per_node;
        item = item * pct.scalapack_nprow * pct.scalapack_npcol;

        int ppp;
        MPI_Comm_size(pct.scalapack_comm, &ppp);

        MPI_Bcast (eigs, num_states, MPI_DOUBLE, item, pct.scalapack_comm);
    }



    delete [] Cij;
    delete [] distCij;
    delete [] distSij;
    delete [] distBij;
    delete [] distAij;
#endif
}
