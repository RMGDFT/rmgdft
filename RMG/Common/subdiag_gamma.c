/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/subdiag_mpi.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void subdiag(STATE *states, rmg_double_t *vh, rmg_double_t *vnuc, rmg_double_t *vxc)
 *   This is the version of subdiag for message passing machines(MPI).
 *   Subspace diagonalizer for the Mehrstellen discretization of the
 *   Kohn-Sham Hamiltonian.
 *
 *    1. Computes <i |A + BV| j> and <i |B| j>
 *    2. Generates the inverse of <i |B| j> and then multiplies <i |A + BV| j>
 *       by this inverse.
 *    3. Diagonalizes the resulting matrix and then remixes the eigenvectors
 *       accordingly.
 * INPUTS
 *   states: points to orbital structure (see main.h)
 *   vh:     Hartree potential
 *   vnuc:   pseudopotential
 *   vxc:    exchange correlation potential
 * OUTPUT
 *   wave functions and eigenvalues in states are updated
 * PARENTS
 *   cdfastrlx.c fastrlx.c moldyn.c init.c quench.c
 * CHILDREN
 *   global_sums.c app_cir.c app_nl.c genvpsi.c pack_ptos.c app_cilr.c 
 *   gather_psi.c app_grad.c
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"
#include "blas.h"


#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

int rmg_dsygvd_gpu(int n, rmg_double_t *a, int lda, rmg_double_t *b, int ldb, 
		   rmg_double_t *w, rmg_double_t *work, int lwork, int *iwork, int liwork, rmg_double_t *wa);

static rmg_double_t *tmp_arrayR;
static rmg_double_t *tmp_array2R;

// Array storage for folded spectrum diagonalization
static int *fs_eigstart;
static int *fs_eigstop;
static int *fs_eigcounts;

/* This subspace diagonalization function uses Scalapack libraries  */

#if GAMMA_PT

/*This function works for gamma point only. The matrix multiplications used to setup matrices
 * A, B and S cannot be used for non-gamma point calculations, that would require
 * change in wavefunctions memory layout - we need to have real and complex parts to form 
 * continous array, but currently our wafunction array consists of complex and real 
 * parts in alternating order. At this time, I do not want to spend time on non-gamma point 
 * calculation, which I do not really need.
 * 
 * I have written non-gamma point versions of the functions subdiag_app_A, subdiag_app_S, 
 * subdiag_app_B, what is required for this function to work in the non-gamma point mode
 * is calculating matrices A, B and S and solve complex generalized eigenvalue
 * problem, i.e use complex version of PDSYGVX (PZHEGVX, I think).*/


#if GAMMA_PT
static void subdiag2_mpi (rmg_double_t * Aij, rmg_double_t * base_mem, rmg_double_t * tmp_psi);
void subdiag_gamma_scalapack (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc);
void subdiag_gamma_lapack (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc);
void subdiag_gamma_magma (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc);
void subdiag_gamma_elpa (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc);
#else
static void subdiag2_mpi (rmg_double_t * Aij, rmg_double_t * base_mem);
#endif
static void symmetrize_matrix (rmg_double_t * matrix, rmg_double_t * unity_matrix, int size, int *desca,
                               int local_size);
static void print_matrix2 (rmg_double_t * matrix, int size);
static void print_dist_matrix (rmg_double_t * dist_matrix, int global_size, int *desca);



// MPI operations are faster on some systems when the memory is allocated by MPI_Alloc_mem
// instead of the normal system malloc.
static rmg_double_t *distAij, *distBij, *distCij, *distIij, *distSij, *distTij;
rmg_double_t *global_matrix;
void init_subdiag(void)
{

    int dist_length, dist_stop, pbasis, num_states, retval, stop, idx;
    int ione = 1, izero = 0;    /* blas constants */
    rmg_double_t t1, t2;
    rmg_double_t time2;

    /*************************** ScaLapack initialization *************************************/

    time2 = my_crtc ();

    pbasis = get_P0_BASIS();
    num_states = ct.num_states;
    dist_stop = pct.scalapack_max_dist_size;
    stop = num_states * num_states;

    if((ct.subdiag_driver == SUBDIAG_LAPACK) ||
       (ct.subdiag_driver == SUBDIAG_LAPACKFS) ||
       (ct.subdiag_driver == SUBDIAG_MAGMA)  ||
       (ct.subdiag_driver == SUBDIAG_MAGMAFS)) dist_stop=stop;


#if !GAMMA_PT
    dist_stop *= 2;
    stop *= 2;
#endif

    /*Get memory for distributed matrices */
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * dist_stop , MPI_INFO_NULL, &distAij);
    if(retval != MPI_SUCCESS) {
        error_handler("Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * dist_stop , MPI_INFO_NULL, &distBij);
    if(retval != MPI_SUCCESS) {
        error_handler("Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * dist_stop , MPI_INFO_NULL, &distSij);
    if(retval != MPI_SUCCESS) {
        error_handler("Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * dist_stop , MPI_INFO_NULL, &distCij);
    if(retval != MPI_SUCCESS) {
        error_handler("Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * dist_stop , MPI_INFO_NULL, &distIij);
    if(retval != MPI_SUCCESS) {
        error_handler("Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * dist_stop * ct.THREADS_PER_NODE, MPI_INFO_NULL, &distTij);
    if(retval != MPI_SUCCESS) {
        error_handler("Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * stop , MPI_INFO_NULL, &global_matrix);
    if(retval != MPI_SUCCESS) {
        error_handler("Error in MPI_Alloc_mem.\n");
    }

    /*Temporary memory that will be used to calculate matrices and to update wavefunctions */
#if GPU_ENABLED
    cudaHostRegister( global_matrix, sizeof(rmg_double_t) * stop, cudaHostRegisterPortable);

    retval = MPI_Alloc_mem(pbasis * ct.num_states * sizeof(rmg_double_t) , MPI_INFO_NULL, &tmp_arrayR);
    if(retval != MPI_SUCCESS) {
        error_handler("Error in MPI_Alloc_mem.\n");
    }
    cudaHostRegister( tmp_arrayR, pbasis * ct.num_states * sizeof(rmg_double_t), cudaHostRegisterPortable);

    retval = MPI_Alloc_mem(pbasis * ct.num_states * sizeof(rmg_double_t) , MPI_INFO_NULL, &tmp_array2R);
    if(retval != MPI_SUCCESS) {
        error_handler("Error in MPI_Alloc_mem.\n");
    }
    cudaHostRegister( tmp_array2R, pbasis * ct.num_states * sizeof(rmg_double_t), cudaHostRegisterPortable);

#else
    my_malloc (tmp_arrayR, pbasis * ct.num_states, rmg_double_t);
    my_malloc (tmp_array2R, pbasis * ct.num_states, rmg_double_t);
#endif


    // Allocate and initialize distribution arrays for folded spectrum
    my_malloc (fs_eigstart, NPES, int);
    my_malloc (fs_eigstop, NPES, int);
    my_malloc (fs_eigcounts, NPES, int);
    for(idx = 0;idx < NPES;idx++) {
        t1 = (rmg_double_t)ct.num_states;
        t1 = t1 / ((rmg_double_t)NPES);
        t2 = t1 * (rmg_double_t)idx;
        fs_eigstart[idx] = (int)rint(t2);
        fs_eigstop[idx] = (int)rint(t1 + t2);
        fs_eigcounts[idx] = fs_eigstop[idx] - fs_eigstart[idx];
        fs_eigstart[idx] *= ct.num_states;
        fs_eigstop[idx] *= ct.num_states;
        fs_eigcounts[idx] *= ct.num_states;
        
    }


    rmg_timings (DIAG_SCALAPACK_INIT, my_crtc () - time2);
    /********************* Scalapack should be initialized ******************************/

}

void subdiag_gamma (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc)
{

#if GPU_ENABLED
      cuCtxSynchronize();
#endif
    switch(ct.subdiag_driver) {
        case SUBDIAG_LAPACK:
        case SUBDIAG_LAPACKFS:
            subdiag_gamma_lapack(states, vh, vnuc, vxc);
            break;
#if GPU_ENABLED
        case SUBDIAG_MAGMA:
        case SUBDIAG_MAGMAFS:
            subdiag_gamma_magma(states, vh, vnuc, vxc);
            break;
#endif
        default:
            subdiag_gamma_scalapack(states, vh, vnuc, vxc);
    }

}

void subdiag_gamma_scalapack (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc)
{
    int idx, st1, st2, ion, nion, ip, pstop, FP0_BASIS;
	int num_states;
    int stop;
    int kidx;
    rmg_double_t *eigs, *work1R, *sintR;
    int ione = 1, izero = 0;    /* blas constants */
    char *uplo = "l", *jobz = "v";
    ION *iptr;
    SPECIES *sp;

    int info = 0;
    rmg_double_t time1, time2, time3;
    rmg_double_t tmp1;
#if !GAMMA_PT
    rmg_double_t *tmp_arrayI;
#endif
//    rmg_double_t *distAij, *distBij, *distCij, *distIij, *distSij;
    rmg_double_t *vtot, *vtot_eig;
    int dist_length, dist_stop, pbasis;

#if GPU_ENABLED
    cublasStatus_t custat;
    cublasOperation_t cu_transT = CUBLAS_OP_T, cu_transN = CUBLAS_OP_N;
#endif
    rmg_double_t alpha1 = 1.0, beta1 = 0.0;

    num_states = ct.num_states;
    pbasis = get_P0_BASIS();
    FP0_BASIS = get_FP0_BASIS();
    stop = num_states * num_states;

#if GPU_ENABLED

//    cublasSetVectorAsync( pbasis * num_states, sizeof( rmg_double_t ), states[0].psiR, ione, ct.gpu_states, ione, ct.cuda_stream );
    cublasSetVector( pbasis * num_states, sizeof( rmg_double_t ), states[0].psiR, ione, ct.gpu_states, ione );

#endif
        
    time1 = my_crtc ();

    if (pct.gridpe == 0)
        printf ("\n SUBSPACE DIAGONALIZATION");

    kidx = states[0].kidx;


    /* Grab some temporary storage */
#  if !GAMMA_PT
    stop = 2 * num_states * num_states;
#  endif

    /*Get memory for global matrices */
    for(idx = 0;idx < stop;idx++) global_matrix[idx] = 0.0;
    my_malloc (vtot_eig, pbasis, rmg_double_t);
    my_malloc (eigs, num_states, rmg_double_t);
    my_malloc (work1R, ct.num_states * 16 , rmg_double_t);


    /*Get vtot on coarse grid */
    my_malloc (vtot, FP0_BASIS, rmg_double_t);
    for (idx = 0; idx < FP0_BASIS; idx++)
        vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx];
    get_vtot_psi (vtot_eig, vtot, FG_NX);

    /*Release memory for vtot, do it already here since vtot is on fine grid */
    my_free (vtot);



    /*Temporary memory that will be used to calculate matrices and to update wavefunctions */
#if !GAMMA_PT
    my_malloc (tmp_arrayI, pbasis * ct.num_states, rmg_double_t);
#endif


    /*************************** ScaLapack initialization *************************************/

    time2 = my_crtc ();

    if (pct.scalapack_pe)
    {

        /*Length of distributed matrices (different on each processor) */
        dist_length =
            NUMROC (&num_states, &pct.desca[4], &pct.scalapack_myrow, &izero,
                    &pct.scalapack_nprow) * NUMROC (&num_states, &pct.desca[4], &pct.scalapack_mycol,
                                                    &izero, &pct.scalapack_npcol);

        /* Every processor for which pct.scalapack_pe should have some data and so dist_length cannot be 0*/
        if (dist_length == 0)
        error_handler(" function NUMROC returned 0, that should not happen");


        /*This holds number of doubles on each PE */
        dist_stop = dist_length;
#if !GAMMA_PT
        dist_stop *= 2;
#endif
        /* Clear distributed matrices */
        for(idx = 0;idx < dist_stop;idx++) {
            distAij[idx] = 0.0; 
        }
        for(idx = 0;idx < dist_stop;idx++) {
            distBij[idx] = 0.0; 
        }
        for(idx = 0;idx < dist_stop;idx++) {
            distSij[idx] = 0.0; 
        }
        for(idx = 0;idx < dist_stop;idx++) {
            distCij[idx] = 0.0; 
        }
        for(idx = 0;idx < dist_stop;idx++) {
            distIij[idx] = 0.0; 
        }

    }
    rmg_timings (DIAG_SCALAPACK_INIT, my_crtc () - time2);
    /********************* Scalapack should be initialized ******************************/



    /*Get matrix Aij */
    {
        char *trans = "t";
        char *trans2 = "n";
        rmg_double_t alpha = 1.0;
        rmg_double_t beta = 0.0;

        time2 = my_crtc ();

        /*Apply AB operator on each wavefunction 
	 tmp_arrayR:  A|psi> + BV|psi> + B|beta>dnm<beta|psi>
	 tmp_array2R:  B|psi> + B|beta>qnm<beta|psi> */

	subdiag_app_AB (states, tmp_arrayR, tmp_array2R, vtot_eig);

        time3 = my_crtc ();
        rmg_timings (DIAG_APP_A, time3 - time2);

#if GPU_ENABLED

        cublasSetVector( pbasis * num_states, sizeof( rmg_double_t ), tmp_arrayR, ione, ct.gpu_temp, ione );

        cublasDgemm(ct.cublas_handle, cu_transT, cu_transN, num_states, num_states, pbasis,
             &alpha, ct.gpu_states, pbasis,
             ct.gpu_temp, pbasis,
             &beta,  ct.gpu_global_matrix, num_states );
        
        cublasGetVector( num_states * num_states, sizeof( rmg_double_t ), ct.gpu_global_matrix, ione, global_matrix, ione );

#else
        /*Global matrix will hold global A matrix */
        dgemm (trans, trans2, &num_states, &num_states, &pbasis, &alpha, states[0].psiR, &pbasis,
               tmp_arrayR, &pbasis, &beta, global_matrix, &num_states);
#endif


        time2 = my_crtc ();
        rmg_timings (DIAG_DGEMM, time2 - time3);

        // Reduce and distribute Aij
        reduce_and_dist_matrix(num_states, global_matrix, distAij, distTij);


        // Now deal with the S operator
        time3 = my_crtc ();
        alpha = ct.vel;
#if GPU_ENABLED
        cublasSetVector( pbasis * num_states, sizeof( rmg_double_t ), pct.ns, ione, ct.gpu_temp, ione );
        cublasDgemm(ct.cublas_handle, cu_transT, cu_transN, num_states, num_states, pbasis,
             &alpha, ct.gpu_states, pbasis,
             ct.gpu_temp, pbasis,
             &beta,  ct.gpu_global_matrix, num_states );

        cublasGetVector( num_states * num_states, sizeof( rmg_double_t ), ct.gpu_global_matrix, ione, global_matrix, ione );

#else
        dgemm (trans, trans2, &num_states, &num_states, &pbasis, &alpha, states[0].psiR, &pbasis,
               pct.ns, &pbasis, &beta, global_matrix, &num_states);
#endif

        time2 = my_crtc ();
        rmg_timings (DIAG_DGEMM, time2 - time3);

        // Reduce and distribute Sij
        reduce_and_dist_matrix(num_states, global_matrix, distSij, distTij);

        /* Apply B operator on each wavefunction */

        time3 = my_crtc ();
        alpha = ct.vel;

		//subdiag_app_B (states, tmp_array2R);
       // for(idx = 0; idx < 100; idx++) printf("\n iaaa  %d  %f ", idx, tmp_array2R[idx]);

#if GPU_ENABLED
        cublasSetVector( pbasis * num_states, sizeof( rmg_double_t ), tmp_array2R, ione, ct.gpu_temp, ione );

        cublasDgemm(ct.cublas_handle, cu_transT, cu_transN, num_states, num_states, pbasis,
             &alpha, ct.gpu_states, pbasis,
             ct.gpu_temp, pbasis,
             &beta,  ct.gpu_global_matrix, num_states );

        cublasGetVector( num_states * num_states, sizeof( rmg_double_t ), ct.gpu_global_matrix, ione, global_matrix, ione );

#else
        dgemm (trans, trans2, &num_states, &num_states, &pbasis, &alpha, states[0].psiR, &pbasis,
               tmp_array2R, &pbasis, &beta, global_matrix, &num_states);
#endif


        time2 = my_crtc ();
        rmg_timings (DIAG_DGEMM, time2 - time3);

        // Reduce and distribute Bij
        reduce_and_dist_matrix(num_states, global_matrix, distBij, distTij);



    }

    /*Create distributed unitary matrix */
    for (st1 = 0; st1 < stop; st1++)
        global_matrix[st1] = 0.0;


    for (st1 = 0; st1 < num_states; st1++)
    {
#if GAMMA_PT
        global_matrix[st1 * num_states + st1] = 1.0;
#else
        global_matrix[2 * (st1 * num_states + st1)] = 1.0;
#endif
    }


    time2 = my_crtc ();
    /*Distribute unitary matrix */
    if (pct.scalapack_pe)
        distribute_mat (pct.desca, global_matrix, distCij, &num_states);

    rmg_timings (DIAG_DISTMAT, my_crtc () - time2);

    time2 = my_crtc ();

    /*Scalapack operations are performed here */
    if (pct.scalapack_pe)
    {

        /*keep an extra copy of distributed unitary matrix */
        QMD_dcopy (dist_stop, distCij, ione, distIij, ione);

        /*Get matrix that is inverse to B */
        {
            int ipiv_size, *ipiv;

            ipiv_size =
                NUMROC (&pct.desca[2], &pct.desca[4], &pct.scalapack_myrow, &pct.desca[6],
                        &pct.scalapack_nprow) + pct.desca[4];
            my_calloc (ipiv, ipiv_size, int);

            /*Inverse of B should be in Cij */
#if GAMMA_PT
            PDGESV (&num_states, &num_states, distBij, &ione, &ione, pct.desca, ipiv, distCij, &ione,
                    &ione, pct.desca, &info);
#else
            PZGESV (&num_states, &num_states, distBij, &ione, &ione, pct.desca, ipiv, distCij, &ione,
                    &ione, pct.desca, &info);
#endif

            if (info)
            {
                printf ("\n PE %d: p{d,z}gesv failed, info is %d", pct.gridpe, info);
                error_handler (" p{d,z}gesv failed");
            }

            my_free (ipiv);
        }

#if 0
        if (pct.gridpe == 0)
            printf ("\n\n Inverse Bij");
        print_dist_matrix (distCij, num_states, pct.desca);
#endif


        /*Multiply inverse of B and and A */
        {
            char *trans = "n";
            rmg_double_t alpha[] = { 1.0, 0.0 };
            rmg_double_t beta[] = { 0.0, 0.0 };

            /*B^-1*A */
#if GAMMA_PT
            PDGEMM (trans, trans, &num_states, &num_states, &num_states, alpha,
                    distCij, &ione, &ione, pct.desca, distAij, &ione, &ione, pct.desca, beta, distBij,
                    &ione, &ione, pct.desca);
#else
            PZGEMM (trans, trans, &num_states, &num_states, &num_states, alpha,
                    distCij, &ione, &ione, pct.desca, distAij, &ione, &ione, pct.desca, beta, distBij,
                    &ione, &ione, pct.desca);
#endif

#if 0
            if (pct.gridpe == 0)
                printf ("\n\n Inverse Bij multiplied by Aij, symmetrized");
            print_dist_matrix (distBij, num_states, pct.desca);
#endif


            /*Multiply the result with Sij, result is in distCij */
#if GAMMA_PT
            PDGEMM (trans, trans, &num_states, &num_states, &num_states, alpha,
                    distSij, &ione, &ione, pct.desca, distBij, &ione, &ione, pct.desca, beta, distCij,
                    &ione, &ione, pct.desca);
#else
            PZGEMM (trans, trans, &num_states, &num_states, &num_states, alpha,
                    distSij, &ione, &ione, pct.desca, distBij, &ione, &ione, pct.desca, beta, distCij,
                    &ione, &ione, pct.desca);
#endif

            /*Copy result into Bij */
            QMD_dcopy (dist_stop, distCij, ione, distBij, ione);
        }


#if 1
        /****************** Find Matrix of Eigenvectors *****************************/
        /* Using lwork=-1, pdsyev should return minimum required size for the work array */
        {
            char *range = "a";
            rmg_double_t vx = 0.0;
            rmg_double_t tol = 0.0;
            int eigs_found, eigvs_found;
            rmg_double_t orfac = 0.0;
            int *iwork, *ifail, *iclustr, lwork;
            rmg_double_t *gap, lwork_tmp, *work2;
            int liwork_tmp, liwork;

            my_malloc (ifail, num_states, int);
            my_malloc (iclustr, 2 * pct.scalapack_nprow * pct.scalapack_npcol, int);
            my_malloc (gap, pct.scalapack_nprow * pct.scalapack_npcol, rmg_double_t);
            lwork = -1;
            liwork = -1;


            PDSYGVX (&ione, jobz, range, uplo, &num_states, distBij, &ione, &ione, pct.desca,
                     distSij, &ione, &ione, pct.desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
                     &eigvs_found, eigs, &orfac, distAij, &ione, &ione, pct.desca, &lwork_tmp, &lwork,
                     &liwork_tmp, &liwork, ifail, iclustr, gap, &info);

            if (info)
            {
                printf ("\n PDSYGVX query failed, info is %d", info);
                error_handler ("PDSYGVX query failed");
            }

            /*set lwork and liwork */
            lwork = (int) lwork_tmp + 1;
            liwork = liwork_tmp;

            my_malloc (work2, lwork, rmg_double_t);
            my_malloc (iwork, liwork, int);

            tol = 1e-15;




            PDSYGVX (&ione, jobz, range, uplo, &num_states, distBij, &ione, &ione, pct.desca,
                     distSij, &ione, &ione, pct.desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
                     &eigvs_found, eigs, &orfac, distAij, &ione, &ione, pct.desca, work2, &lwork, iwork,
                     &liwork, ifail, iclustr, gap, &info);

            if (info)
            {
                printf ("\n PDSYGVX failed, info is %d", info);
                error_handler ("PDSYGVX failed");
            }


            /*If subspace diagonalization is used everystep, use eigenvalues obtained here 
             * as the correct eigenvalues*/
            if (ct.diag == 1)
                for (st1 = 0; st1 < ct.num_states; st1++)
                    states[st1].eig[0] = eigs[st1];


            my_free (ifail);
            my_free (iclustr);
            my_free (gap);
            my_free (work2);
            my_free (iwork);

        }
#endif


        /*Gather result onto global_matrix */
        matgather (distAij, pct.desca, global_matrix, num_states);


    }                           /*end if (pct.scalapack_pe) */

    /*Other processors should have global_matrix set to 0 */
    else
        for (idx = 0; idx < stop; idx++)
            global_matrix[idx] = 0.0;


    rmg_timings (DIAG_MATRIX_TIME, (my_crtc () - time2));




    /*Finally, sum global_matrix over all PEs */
    time3 = my_crtc ();

//    global_sums (global_matrix, &stop, pct.grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, global_matrix, stop, MPI_DOUBLE, MPI_SUM, pct.scalapack_comm);
    rmg_timings (DIAG_GLOB_SUMS, my_crtc () - time3);



    /*If some processors did not participate in Scalapack,
     * broadcast eigenvalues, since only Scalapack processors have updated eigenvalues*/
    if ((pct.scalapack_nprow * pct.scalapack_npcol != pct.scalapack_npes) && (ct.diag == 1))
    {

        time2 = my_crtc ();
	int item;
	item = pct.thisimg % ct.images_per_node;
	item = item * pct.scalapack_nprow * pct.scalapack_npcol;

	int ppp;
	MPI_Comm_size(pct.scalapack_comm, &ppp);

	MPI_Bcast (eigs, num_states, MPI_DOUBLE, item, pct.scalapack_comm);

	/*Assign eigenvalues */
	for (st1 = 0; st1 < ct.num_states; st1++)
		states[st1].eig[0] = eigs[st1];

	rmg_timings (DIAG_BCAST_EIGS, (my_crtc () - time2));
    }




    time2 = my_crtc ();

    /* Do the orbital update in here */
#if GAMMA_PT
#if GPU_ENABLED

    cublasSetVector( num_states * num_states, sizeof( rmg_double_t ), global_matrix, ione, ct.gpu_global_matrix, ione );

    alpha1 = 1.0;
    beta1 = 0.0;
    custat = cublasDgemm(ct.cublas_handle, cu_transN, cu_transN, pbasis, num_states, num_states,
		    &alpha1, 
		    ct.gpu_states, pbasis,
		    ct.gpu_global_matrix, num_states,
		    &beta1,  ct.gpu_temp, pbasis );

    cublasGetVector( pbasis * num_states, sizeof( rmg_double_t ), ct.gpu_temp, ione, states->psiR, ione );

#else
    subdiag2_mpi (global_matrix, states->psiR, tmp_arrayR);
#endif
#else
    subdiag2_mpi (global_matrix, states->psiR);
#endif


    rmg_timings (DIAG_WAVEUP_TIME, (my_crtc () - time2));


    /* release our temporary storage */

#if !GAMMA_PT
    my_free (tmp_arrayI);
#endif
    my_free (work1R);
    my_free (eigs);
    my_free (vtot_eig);

    rmg_timings (DIAG_TIME, (my_crtc () - time1));

}




#  if GAMMA_PT
/* This routine is used to do the subspace rotation of the orbitals. Each
 * thread handles a specific portion of the real space domain.
 */
static void subdiag2_mpi (rmg_double_t * Aij, rmg_double_t * base_mem, rmg_double_t * tmp_psi)
{
	int idx;

	char *trans = "n";
	rmg_double_t alpha = 1.0;
	rmg_double_t beta = 0.0;
	int pbasis = get_P0_BASIS();
	int num_states = ct.num_states;

	dgemm (trans, trans, &pbasis, &num_states, &num_states, &alpha, base_mem, &pbasis, Aij,
			&num_states, &beta, tmp_psi, &pbasis);

	for (idx = 0; idx < num_states * pbasis; idx++)
		base_mem[idx] = tmp_psi[idx];


}                               /* end subdiag2_mpi */

#  else

/* This routine is used to do the subspace rotation of the orbitals. Each
 * thread handles a specific portion of the real space domain.
 */
void subdiag2_mpi (rmg_double_t * Aij, rmg_double_t * base_mem)
{
	int idx, st1, st2, P0_BASIS;
	rmg_double_t *rptr;
	rmg_double_t *work1R, *work2R;
	rmg_double_t *work1I, *work2I;

	my_malloc (work1R, 4 * ct.num_states, rmg_double_t);
	work2R = work1R + ct.num_states;
	work1I = work2R + ct.num_states;
	work2I = work1I + ct.num_states;

	rptr = base_mem;

        P0_BASIS = get_P0_BASIS();
	for (idx = 0; idx < P0_BASIS(); idx++)
	{

		/* We make a temporary copy and store it in work2 otherwise the
		 * cache thrashing is horrible on the O2K.
		 */
		for (st2 = 0; st2 < ct.num_states; st2++)
		{
			work2R[st2] = rptr[2 * st2 * P0_BASIS + idx];
			work2I[st2] = rptr[2 * st2 * P0_BASIS + P0_BASIS + idx];
		}

		for (st1 = 0; st1 < ct.num_states; st1++)
		{

			work1R[st1] = 0.0;
			work1I[st1] = 0.0;
			for (st2 = 0; st2 < ct.num_states; st2++)
			{

				work1R[st1] += Aij[2 * (st1 * ct.num_states + st2)] * work2R[st2];
				work1R[st1] -= Aij[2 * (st1 * ct.num_states + st2) + 1] * work2I[st2];
				work1I[st1] += Aij[2 * (st1 * ct.num_states + st2)] * work2I[st2];
				work1I[st1] += Aij[2 * (st1 * ct.num_states + st2) + 1] * work2R[st2];

			}                   /* st2 */

		}                       /* end for */

		/* update all wavefunctions for this *idx* */
		for (st1 = 0; st1 < ct.num_states; st1++)
		{
			rptr[2 * st1 * P0_BASIS + idx] = work1R[st1];
			rptr[2 * st1 * P0_BASIS + P0_BASIS + idx] = work1I[st1];
		}

	}                           /* idx */

	my_free (work1R);

}                               /* end subdiag2_mpi */


#  endif

static void print_matrix2 (rmg_double_t * matrix, int size)
{
	int i, j;

	//printf("\n\n");

	for (i = 0; i < size; i++)
	{
		printf ("\n");
#if GAMMA_PT
		for (j = 0; j < size; j++)
			printf ("%g  ", matrix[i * size + j]);
#else
		for (j = 0; j < 2 * size; j += 2)
			printf ("%  f %f  ", matrix[2 * i * size + j], matrix[2 * i * size + j + 1]);
#endif
	}

}

static void print_dist_matrix (rmg_double_t * dist_matrix, int global_size, int *desca)
{
	rmg_double_t *glob_matrix;
	int stop;

	stop = global_size * global_size;

	my_calloc (glob_matrix, stop, rmg_double_t);


	if (pct.scalapack_pe)
		matgather (dist_matrix, desca, glob_matrix, global_size);


	/*Sum Aij over all PEs */
	global_sums (glob_matrix, &stop, pct.grid_comm);

	if (!pct.gridpe)
		print_matrix2 (glob_matrix, global_size);

	my_free (glob_matrix);


}



/*This works with distributed matrices*/
static void symmetrize_matrix (rmg_double_t * matrix, rmg_double_t * unity_matrix, int size, int *desca,
		int local_size)
{
	int stop, ione = 1;
	rmg_double_t *temp_unity_matrix, *temp_matrix;
	rmg_double_t alpha[] = { 0.5, 0.0 };
	rmg_double_t beta[] = { 0.5, 0.0 };
	char *trans = "n";
#if GAMMA_PT
	char *trans2 = "t";
#else
	char *trans2 = "c";
#endif

	stop = local_size;
#if !GAMMA_PT
	stop *= 2;
#endif


	/*Get memory */
	my_malloc (temp_matrix, stop, rmg_double_t);
	my_calloc (temp_unity_matrix, stop, rmg_double_t);

	/*Copy matrix into temp_matrix */
	QMD_dcopy (stop, matrix, ione, temp_matrix, ione);

	/*Local copy of unity matrix, this is done so that the unitary matrix that was passed here does not change */
	QMD_dcopy (stop, unity_matrix, ione, temp_unity_matrix, ione);


	/*Symmetric (or Hermitian) matrix will be obtained as
	 * A = 0.5(I*A + A^T)*/

#if GAMMA_PT
	PDGEMM (trans, trans2, &size, &size, &size, alpha,
			temp_unity_matrix, &ione, &ione, desca, temp_matrix, &ione, &ione, desca, beta, matrix,
			&ione, &ione, desca);
#else
	PZGEMM (trans, trans2, &size, &size, &size, alpha,
			temp_unity_matrix, &ione, &ione, desca, temp_matrix, &ione, &ione, desca, beta, matrix,
			&ione, &ione, desca);
#endif


	/*Release memory */
	my_free (temp_matrix);
	my_free (temp_unity_matrix);

}

#endif



#if GAMMA_PT
void subdiag_gamma_lapack (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc)
{
	int idx, st1, st2, ion, nion, ip, pstop, P0_BASIS, FP0_BASIS;
	int num_states;
	int stop;
	int kidx;
	rmg_double_t *eigs, *work1R, *sintR;
	int ione = 1, izero = 0;    /* blas constants */
	char *uplo = "l", *jobz = "v";
	ION *iptr;
	SPECIES *sp;

	int info = 0;
	rmg_double_t time1, time2, time3;
	rmg_double_t tmp1;
	rmg_double_t *vtot, *vtot_eig;
	int dist_stop, pbasis;

#if GPU_ENABLED
	cublasStatus_t custat;
	cublasOperation_t cu_transT = CUBLAS_OP_T, cu_transN = CUBLAS_OP_N;
#endif
	rmg_double_t alpha1 = 1.0, beta1 = 0.0;

	num_states = ct.num_states;
	pbasis = get_P0_BASIS();
        P0_BASIS = pbasis;
        FP0_BASIS = get_FP0_BASIS();

	stop = num_states * num_states;

#if GPU_ENABLED

	cublasSetVector( pbasis * num_states, sizeof( rmg_double_t ), states[0].psiR, ione, ct.gpu_states, ione );

#endif

	time1 = my_crtc ();

	if (pct.gridpe == 0)
		printf ("\n SUBSPACE DIAGONALIZATION");

	kidx = states[0].kidx;


	/*Get memory for global matrices */
	for(idx = 0;idx < stop;idx++) global_matrix[idx] = 0.0;
	my_malloc (vtot_eig, P0_BASIS, rmg_double_t);
	my_malloc (eigs, num_states, rmg_double_t);
	my_malloc (work1R, ct.num_states * 16 , rmg_double_t);


	/*Get vtot on coarse grid */
	my_malloc (vtot, FP0_BASIS, rmg_double_t);
	for (idx = 0; idx < FP0_BASIS; idx++)
		vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx];
	get_vtot_psi (vtot_eig, vtot, FG_NX);

	/*Release memory for vtot, do it already here since vtot is on fine grid */
	my_free (vtot);





	/*************************** ScaLapack initialization *************************************/

	time2 = my_crtc ();

	/*This holds number of doubles on each PE */
	dist_stop = stop;
	/* Clear distributed matrices */
	for(idx = 0;idx < dist_stop;idx++) {
		distAij[idx] = 0.0; 
	}
	for(idx = 0;idx < dist_stop;idx++) {
		distBij[idx] = 0.0; 
	}
	for(idx = 0;idx < dist_stop;idx++) {
		distSij[idx] = 0.0; 
	}
	for(idx = 0;idx < dist_stop;idx++) {
		distCij[idx] = 0.0; 
	}

	rmg_timings (DIAG_SCALAPACK_INIT, my_crtc () - time2);
	/********************* Scalapack should be initialized ******************************/



	/*Get matrix Aij */
	{
		char *trans = "t";
		char *trans2 = "n";
		rmg_double_t alpha = 1.0;
		rmg_double_t beta = 0.0;

		time2 = my_crtc ();

		/*Apply A operator on each wavefunction 
		 * S operator is also applied, the result is returned in tmp_array2R*/
		subdiag_app_AB (states, tmp_arrayR, tmp_array2R, vtot_eig);

		time3 = my_crtc ();
		rmg_timings (DIAG_APP_A, time3 - time2);

#if GPU_ENABLED

		cublasSetVector( pbasis * num_states, sizeof( rmg_double_t ), tmp_arrayR, ione, ct.gpu_temp, ione );

		cublasDgemm(ct.cublas_handle, cu_transT, cu_transN, num_states, num_states, pbasis,
				&alpha, ct.gpu_states, pbasis,
				ct.gpu_temp, pbasis,
				&beta,  ct.gpu_global_matrix, num_states );

		cublasGetVector( num_states * num_states, sizeof( rmg_double_t ), ct.gpu_global_matrix, ione, global_matrix, ione );

#else
		/*Global matrix will hold global A matrix */
		dgemm (trans, trans2, &num_states, &num_states, &pbasis, &alpha, states[0].psiR, &pbasis,
				tmp_arrayR, &pbasis, &beta, global_matrix, &num_states);
#endif


		time2 = my_crtc ();
		rmg_timings (DIAG_DGEMM, time2 - time3);

		// Reduce Aij
		//global_sums (global_matrix, &stop, pct.grid_comm);
		MPI_Allreduce(MPI_IN_PLACE, global_matrix, stop, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
		rmg_timings (DIAG_GLOB_SUMS, my_crtc () - time2);
		QMD_dcopy (stop, global_matrix, ione, distAij, ione);


		// Now deal with the S operator
		time3 = my_crtc ();
		alpha = ct.vel;
#if GPU_ENABLED
		cublasSetVector( pbasis * num_states, sizeof( rmg_double_t ), pct.ns, ione, ct.gpu_temp, ione );
		cublasDgemm(ct.cublas_handle, cu_transT, cu_transN, num_states, num_states, pbasis,
				&alpha, ct.gpu_states, pbasis,
				ct.gpu_temp, pbasis,
				&beta,  ct.gpu_global_matrix, num_states );

		cublasGetVector( num_states * num_states, sizeof( rmg_double_t ), ct.gpu_global_matrix, ione, global_matrix, ione );

#else
		dgemm (trans, trans2, &num_states, &num_states, &pbasis, &alpha, states[0].psiR, &pbasis,
				pct.ns, &pbasis, &beta, global_matrix, &num_states);
#endif

		time2 = my_crtc ();
		rmg_timings (DIAG_DGEMM, time2 - time3);

		// Reduce Sij
		//global_sums (global_matrix, &stop, pct.grid_comm);
		MPI_Allreduce(MPI_IN_PLACE, global_matrix, stop, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
		rmg_timings (DIAG_GLOB_SUMS, my_crtc () - time2);
		QMD_dcopy (stop, global_matrix, ione, distSij, ione);

		/* Apply B operator on each wavefunction */
		time2 = my_crtc ();
//		subdiag_app_B (states, tmp_array2R);

		time3 = my_crtc ();
		alpha = ct.vel;
		rmg_timings (DIAG_APP_B, time3 - time2);

#if GPU_ENABLED
		cublasSetVector( pbasis * num_states, sizeof( rmg_double_t ), tmp_array2R, ione, ct.gpu_temp, ione );

		cublasDgemm(ct.cublas_handle, cu_transT, cu_transN, num_states, num_states, pbasis,
				&alpha, ct.gpu_states, pbasis,
				ct.gpu_temp, pbasis,
				&beta,  ct.gpu_global_matrix, num_states );

		cublasGetVector( num_states * num_states, sizeof( rmg_double_t ), ct.gpu_global_matrix, ione, global_matrix, ione );

#else
		dgemm (trans, trans2, &num_states, &num_states, &pbasis, &alpha, states[0].psiR, &pbasis,
				tmp_array2R, &pbasis, &beta, global_matrix, &num_states);
#endif


		time2 = my_crtc ();
		rmg_timings (DIAG_DGEMM, time2 - time3);

		// Reduce Bij
		//global_sums (global_matrix, &stop, pct.grid_comm);
		MPI_Allreduce(MPI_IN_PLACE, global_matrix, stop, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
		rmg_timings (DIAG_GLOB_SUMS, my_crtc () - time2);
		QMD_dcopy (stop, global_matrix, ione, distBij, ione);


	}

	/*Create unitary matrix */
	for (st1 = 0; st1 < stop; st1++)
		distCij[st1] = 0.0;


	for (st1 = 0; st1 < num_states; st1++)
	{
		distCij[st1 * num_states + st1] = 1.0;
	}

	time2 = my_crtc ();


	/*Get matrix that is inverse to B */
	{
		int *ipiv;

		my_calloc (ipiv, num_states, int);

		/*Inverse of B should be in Cij */
		dgesv (&num_states, &num_states, global_matrix, &num_states, ipiv, distCij, &num_states, &info);

		if (info)
		{
			printf ("\n PE %d: p{d,z}gesv failed, info is %d", pct.gridpe, info);
			error_handler (" p{d,z}gesv failed");
		}

		my_free (ipiv);
	}

	/*Multiply inverse of B and and A */
	{
		char *trans = "n";

		/*B^-1*A */
		dgemm (trans, trans, &num_states, &num_states, &num_states, &alpha1,
				distCij, &num_states, distAij, &num_states, &beta1, distBij,
				&num_states);


		/*Multiply the result with Sij, result is in distCij */
		dgemm (trans, trans, &num_states, &num_states, &num_states, &alpha1,
				distSij, &num_states, distBij, &num_states, &beta1, distCij,
				&num_states);

	}


#if 1
	/****************** Find Matrix of Eigenvectors *****************************/
	/* Using lwork=-1, pdsyev should return minimum required size for the work array */
	{
		char *range = "A";
		rmg_double_t vx = 0.0;
		rmg_double_t tol = 0.0;
		int eigs_found, eigvs_found;
		rmg_double_t orfac = 0.0;
		int *iwork, *ifail, lwork;
		rmg_double_t *gap, lwork_tmp, *work2;
		int liwork_tmp, liwork;

		my_malloc (ifail, num_states, int);
		lwork = -1;
		liwork = -1;

		lwork = 2 * num_states * num_states + 6 * num_states + 2;
		liwork = 6*num_states;
		my_malloc (work2, lwork, rmg_double_t);
		my_malloc (iwork, liwork, int);

		tol = 1e-15;

                if(ct.subdiag_driver == SUBDIAG_LAPACKFS) {

                    rmg_double_t *hwork;
                    my_malloc(hwork, lwork, rmg_double_t);
                    info = rmg_folded_spectrum_cpu(num_states, distCij, num_states, distSij, num_states,
                                eigs, hwork, lwork, iwork, liwork, distAij);
                    my_free(hwork);
                    QMD_dcopy(num_states*num_states, distCij, 1, global_matrix, 1);

                }
                else {

                    dsygvx_  (&ione, jobz, range, uplo, &num_states, distCij, &num_states, distSij, &num_states,
                        &vx, &vx, &ione, &ione,  &tol, &eigs_found, eigs, global_matrix, &num_states, work2, 
                        &lwork, iwork, ifail, &info);

                }

		if (info)
		{
			printf ("\n PDSYGVX failed, info is %d", info);
			error_handler ("PDSYGVX failed");
		}


		/*If subspace diagonalization is used everystep, use eigenvalues obtained here 
		 * as the correct eigenvalues*/
		if (ct.diag == 1)
			for (st1 = 0; st1 < ct.num_states; st1++)
				states[st1].eig[0] = eigs[st1];


		my_free (ifail);
		my_free (work2);
		my_free (iwork);

	}
#endif


	rmg_timings (DIAG_MATRIX_TIME, (my_crtc () - time2));




	/*Finally, sum global_matrix over all PEs */
	time3 = my_crtc ();

	rmg_timings (DIAG_GLOB_SUMS, my_crtc () - time3);


	time2 = my_crtc ();

	/* Do the orbital update in here */
#if GPU_ENABLED

	cublasSetVector( num_states * num_states, sizeof( rmg_double_t ), global_matrix, ione, ct.gpu_global_matrix, ione );

	alpha1 = 1.0;
	beta1 = 0.0;
	custat = cublasDgemm(ct.cublas_handle, cu_transN, cu_transN, pbasis, num_states, num_states,
			&alpha1, 
			ct.gpu_states, pbasis,
			ct.gpu_global_matrix, num_states,
			&beta1,  ct.gpu_temp, pbasis );

	cublasGetVector( pbasis * num_states, sizeof( rmg_double_t ), ct.gpu_temp, ione, states->psiR, ione );

#else
	subdiag2_mpi (global_matrix, states->psiR, tmp_arrayR);
#endif


	rmg_timings (DIAG_WAVEUP_TIME, (my_crtc () - time2));


	/* release our temporary storage */

	my_free (work1R);
	my_free (eigs);
	my_free (vtot_eig);

	rmg_timings (DIAG_TIME, (my_crtc () - time1));

} // end subdiag_gamma_lapack



int rmg_folded_spectrum_cpu(int n, rmg_double_t *a, int lda, rmg_double_t *b, int ldb, 
		rmg_double_t *w, rmg_double_t *work, int lwork, int *iwork, int liwork, rmg_double_t *wa)
{

    int ione=1, itype=1, info=0, idx;
    rmg_double_t rone = 1.0;
    int ix, iy, eig_index, its, n_win, n_start, st, st1, i_width, ireps;
    char *trans="n", *transt="t", *transn="n";
    char *cuplo = "l", *side="l", *diag="n", *jobz="V";
    rmg_double_t alpha, beta=0.0, lambda, *d_p0, *d_p1, t1, t2;
    rmg_double_t *Vdiag, sum, *V, *G, *p0, *p1, *n_eigs, *tarr, *darr, *sarr, *Asave;
    rmg_double_t time1, time2, time3, r_width;
    int eig_step, eig_start, eig_stop, map, istride, ibase, omp_tid;
     
    // Folded spectrum method is parallelized over PE's. Each PE gets assigned
    // a subset of the eigenvectors.
    t1 = (rmg_double_t)n;
    t1 = t1 / ((rmg_double_t)NPES);
    t2 = t1 * (rmg_double_t)pct.gridpe;
    eig_start = (int)rint(t2);
    eig_stop = (int)rint(t1 + t2);
    eig_step = eig_stop - eig_start;
    if(pct.gridpe == (NPES - 1)) eig_stop = n;

    // Set width of window in terms of a percentage of n. Larger values will be slower but
    // exhibit behavior closer to full diagonalization.
    r_width = 0.3;
    t1 = (rmg_double_t)n;
    n_win = (int)(r_width * t1);

    // Find start of interval
    ix = n_win - eig_step;
    if(ix < 4)
        error_handler("Too few PE's to use folded spectrum method for this problem");
    if(ix % 2) {
        n_win++;
        ix = n_win - eig_step;
    }

    n_start = eig_start - ix/2;
    if(n_start < 0)n_start = 0;
    if((n_start + n_win) > n) {
        n_start = n - n_win;
    }

    my_malloc(Vdiag, n, rmg_double_t);
    my_malloc(p0, n, rmg_double_t);
    my_malloc(p1, n, rmg_double_t);
    n_eigs = distTij;

    my_malloc(d_p0, n, rmg_double_t);
    my_malloc(d_p1, n, rmg_double_t);
    my_malloc(Asave, n*n, rmg_double_t);
    my_malloc(tarr, n, rmg_double_t);

    time1=my_crtc();
    //  Form a Cholesky factorization of B.
    dpotrf(cuplo, &n, b, &ldb, &info);
    if( info != 0 ) {
        error_handler("dpotrf failure");
    }

    time2=my_crtc();
    Dprintf("DPOTRF1  = %12.6f",time2-time1);
    time1=my_crtc();

    //  Transform problem to standard eigenvalue problem
    dsygst_(&itype, cuplo, &n, a, &lda, b, &ldb, &info);
    if( info != 0 ) {
        error_handler("dsygst failure");
    }

    time2=my_crtc();
    Dprintf("DSYGST  = %12.6f",time2-time1);
    time1=my_crtc();


    // We need to wait until a is diagonally dominant so we skip the first 3 steps
    if(ct.scf_steps < 3) {

        dsyevd_(jobz, cuplo, &n, a, &lda, w,
                         work,  &lwork,
                         iwork, &liwork,
                         &info);
        if( info != 0 ) {
            error_handler("dsyevd failure");
        }
        //   For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
        //        backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
        //   dtrsm_( "Leftx", uplo, trans, "Non-unit", &n, &n, &rone, b, &ldb, a, &lda );
        //
        dtrsm_ (side, cuplo, transt, diag, &n, &n, &rone, b, &ldb, a, &lda);

    }
    else {
       
        V=tmp_arrayR;
        G=tmp_array2R;

        // AX=lambdaX  store a copy of A in Asave
        QMD_dcopy (n*n, a, 1, Asave, 1);

        // Zero out matrix of eigenvectors and eigenvalues
        for(idx=0;idx < n*n;idx++) V[idx] = 0.0;
        for(idx=0;idx < n;idx++) n_eigs[idx] = 0.0;
     

        time1=my_crtc();

        // Do the submatrix along the diagonal to get starting values for folded spectrum
        //--------------------------------------------------------------------
        for(ix = 0;ix < n_win;ix++){
            for(iy = 0;iy < n_win;iy++){
                G[ix*n_win + iy] = Asave[(n_start+ix)*n + n_start + iy];
            }
        }
        QMD_dcopy (n_win * n_win, G, 1, a, 1);
        dsyevd_(jobz, cuplo, &n_win, a, &n_win, &w[n_start],
                        work, &lwork,
                        iwork, &liwork,
                        &info);
        if( info != 0 ) {
                error_handler("dsyevd failure");
        }
        //--------------------------------------------------------------------

        QMD_dcopy (n_win * n_win, a, 1, G, 1);

        for(ix = 0;ix < n_win;ix++) {
            Vdiag[ix] = 1.0;
            if(G[ix*n_win + ix] < 0.0) Vdiag[ix] = -1.0;
        }

        // Store the eigen vector from the submatrix
        for(ix=0;ix<n_win;ix++) {
            if(((n_start+ix) >= eig_start) && ((n_start+ix) < eig_stop)) {
                for(iy=0;iy<n_win;iy++) {
                      V[(ix + n_start)*n + n_start + iy] = Vdiag[ix] * G[ix * n_win + iy];
                }
            }
        }


        time2=my_crtc();
        dprintf("SUBMATRIX = %12.6f",time2-time1);
        time1=my_crtc();

        // Apply folded spectrum to this PE's range of eigenvectors
        for(eig_index = eig_start;eig_index < eig_stop;eig_index++) {

                lambda = w[eig_index];
                n_eigs[eig_index] = lambda;

                for(ix=0;ix<ct.num_states;ix++){
                    Asave[ix*ct.num_states + ix] -= lambda;
                }
                alpha = 1.0;
                QMD_dcopy (ct.num_states*ct.num_states, Asave, 1, global_matrix, 1);

                // Restore matrix for next pass
                for(ix=0;ix<ct.num_states;ix++){
                    Asave[ix*ct.num_states + ix] += lambda;
                }

                QMD_dcopy (n, &V[eig_index*n], 1, d_p0, 1);

                alpha = -0.5;
                beta = 0.0;
                for(its = 0;its < 6;its++) {
                    dgemv_(transn, &n, &n, &alpha, global_matrix, &n, d_p0, &ione, &beta, d_p1, &ione);
                    daxpy_(&n, &rone, d_p1, &ione, d_p0, &ione);
                }
                // Renormalize
                t1 = dnrm2_(&n, d_p0, &ione);
                t1 = 1.0 / t1;
                dscal_(&n, &t1, d_p0, &ione);

                QMD_dcopy (n, d_p0, 1, &V[eig_index*n], 1);
        }

        time2=my_crtc();
        dprintf("FOLDED SPECTRUM = %12.6f",time2-time1);
        time1=my_crtc();

        // Make sure all PE's have all eigenvectors. Possible optimization here would be to 
        // overlap computation in the above loop with communication here.
//        MPI_Allreduce(MPI_IN_PLACE, V, n*n, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        MPI_Allgatherv(MPI_IN_PLACE, eig_step*n, MPI_DOUBLE, V, fs_eigcounts, fs_eigstart, MPI_DOUBLE, pct.grid_comm);

        time2=my_crtc();
        dprintf("MPI_ALLREDUCE1  = %12.6f",time2-time1);
        time1=my_crtc();

        // Do the same for the eigenvalues
        // Could replace this with an MPI_Allgatherv as well
        MPI_Allreduce(MPI_IN_PLACE, n_eigs, n, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

        time2=my_crtc();
        dprintf("MPI_ALLREDUCE2  = %12.6f",time2-time1);

        // Copy summed eigs back to w
        for(idx = 0;idx < n;idx++) w[idx] = n_eigs[idx];

        time1=my_crtc();

        // Gram-Schmidt ortho for eigenvectors.
        alpha = 1.0;
        beta = 0.0;

        // Overlaps
        QMD_dcopy (n*n, V, 1, a, 1);
        dsyrk_ (cuplo, transt, &n, &n, &alpha, a, &n, &beta, global_matrix, &n);
        time2=my_crtc();
        dprintf("OVERLAPS  = %12.6f",time2-time1);
        time1=my_crtc();

        // Cholesky factorization
        dpotrf(cuplo, &n, global_matrix, &n, &info);

        time2=my_crtc();
        dprintf("CHOLESKY  = %12.6f",time2-time1);
        time1=my_crtc();

        // Get inverse of diagonal elements
        for(ix = 0;ix < n;ix++) tarr[ix] = 1.0 / global_matrix[n * ix + ix];

//----------------------------------------------------------------
        for(idx = 0;idx < n*n;idx++)G[idx] = 0.0;
#pragma omp parallel private(idx,st,st1,omp_tid,sarr)
{
        omp_tid = omp_get_thread_num();
        if(omp_tid == 0) my_malloc(darr, n * omp_get_num_threads(), rmg_double_t);
#pragma omp barrier

#pragma omp for schedule(static, 1) nowait
        for(idx = eig_start;idx < eig_stop;idx++) {

            sarr = &darr[omp_tid*n];

            for (st = 0; st < n; st++) sarr[st] = V[st*n + idx];

            for (st = 0; st < n; st++) {

                sarr[st] *= tarr[st];

                for (st1 = st+1; st1 < n; st1++) {
                    sarr[st1] -= global_matrix[st1 + n*st] * sarr[st];
                }

            }

            for (st = 0; st < n; st++) G[st*n + idx] = sarr[st];

        }
} // end omp section

        my_free(darr);

        time2=my_crtc();
        dprintf("GRAM  = %12.6f",time2-time1);
        time1=my_crtc();


        // A matrix transpose here would let us use an Allgatherv which would be
        // almost twice as fast for the communications part.
        MPI_Allreduce(MPI_IN_PLACE, G, n*n, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        QMD_dcopy (n*n, G, 1, a, 1);

        time2=my_crtc();
        dprintf("MPI_ALLREDUCE3  = %12.6f",time2-time1);
        time1=my_crtc();


        dtrsm_ (side, cuplo, transt, diag, &n, &n, &rone, b, &ldb, a, &lda);
        time2=my_crtc();
        dprintf("DTRSM  = %12.6f",time2-time1);

    }

    my_free(tarr);
    my_free(Asave);
    my_free(d_p1);
    my_free(d_p0);
    my_free(p1);
    my_free(p0);
    my_free(Vdiag);
    return 0;
} 
#endif


#if GAMMA_PT
#if GPU_ENABLED
#if MAGMA_LIBS
void subdiag_gamma_magma (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc)
{
	int idx, st1, st2, ion, nion, ip, pstop;
	int num_states;
	int stop, incy;
	int kidx;
	rmg_double_t *eigs, *work1R, *sintR, rzero=0.0, rone=1.0;
	int ione = 1, izero = 0;    /* blas constants */
	char *uplo = "l", *jobz = "V";
	ION *iptr;
	SPECIES *sp;

	int info = 0;
	rmg_double_t time1, time2, time3;
	rmg_double_t tmp1;
	rmg_double_t *vtot, *vtot_eig;
	rmg_double_t *gpuAij, *gpuBij, *gpuCij, *gpuSij;
	int dist_stop, pbasis;
	cublasStatus_t custat;
	cublasOperation_t cu_transT = CUBLAS_OP_T, cu_transN = CUBLAS_OP_N;
	rmg_double_t alpha1 = 1.0, beta1 = 0.0;

	num_states = ct.num_states;
	pbasis =get_P0_BASIS();
	stop = num_states * num_states;


	// Start wavefunctions transferring to the GPU
	cublasSetVector(pbasis * num_states, sizeof( rmg_double_t ), states[0].psiR, ione, ct.gpu_states, ione );

	gpuAij = ct.gpu_work1;
	gpuBij = ct.gpu_work2;
	gpuCij = ct.gpu_work3;
	gpuSij = ct.gpu_work4;

	time1 = my_crtc ();

	if (pct.gridpe == 0)
		printf ("\n SUBSPACE DIAGONALIZATION");

	kidx = states[0].kidx;


	/*Get memory for global matrices */
	for(idx = 0;idx < stop;idx++) global_matrix[idx] = 0.0;
	my_malloc (vtot_eig,get_P0_BASIS(), rmg_double_t);
	my_malloc (eigs, 2*num_states, rmg_double_t);
	my_malloc (work1R, ct.num_states * 16 , rmg_double_t);


	/*Get vtot on coarse grid */
	my_malloc (vtot, get_FP0_BASIS(), rmg_double_t);
	for (idx = 0; idx < get_FP0_BASIS(); idx++)
		vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx];
	get_vtot_psi (vtot_eig, vtot, FG_NX);

	/*Release memory for vtot, do it already here since vtot is on fine grid */
	my_free (vtot);




	/*************************** ScaLapack initialization *************************************/

	time2 = my_crtc ();

	/*This holds number of doubles on each PE */
	dist_stop = stop;
	/* Clear distributed matrices */
	//    for(idx = 0;idx < dist_stop;idx++) {
	//        distSij[idx] = 0.0; 
	//    }
	//    for(idx = 0;idx < dist_stop;idx++) {
	//        distCij[idx] = 0.0; 
	//    }

	rmg_timings (DIAG_SCALAPACK_INIT, my_crtc () - time2);
	/********************* Scalapack should be initialized ******************************/



	/*Get matrix Aij */
	{
		char *trans = "t";
		char *trans2 = "n";
		rmg_double_t alpha = 1.0;
		rmg_double_t beta = 0.0;

		time2 = my_crtc ();

                /*Apply AB operator on each wavefunction 
                 tmp_arrayR:  A|psi> + BV|psi> + B|beta>dnm<beta|psi>
                 tmp_array2R:  B|psi> + B|beta>qnm<beta|psi> */

                subdiag_app_AB (states, tmp_arrayR, tmp_array2R, vtot_eig);

		time3 = my_crtc ();
		rmg_timings (DIAG_APP_A, time3 - time2);


		/*Global matrix will hold global A matrix */
		cublasSetVector(pbasis * num_states, sizeof( rmg_double_t ), tmp_arrayR, ione, ct.gpu_temp, ione );
		cublasDgemm(ct.cublas_handle, cu_transT, cu_transN, num_states, num_states, pbasis,
				&alpha, ct.gpu_states, pbasis,
				ct.gpu_temp, pbasis,
				&beta,  ct.gpu_global_matrix, num_states );
		cublasGetVector(num_states * num_states, sizeof( rmg_double_t ), ct.gpu_global_matrix, ione, global_matrix, ione );

		time2 = my_crtc ();
		rmg_timings (DIAG_DGEMM, time2 - time3);

		// Reduce Aij and store it on the GPU
		//global_sums (global_matrix, &stop, pct.grid_comm);
		MPI_Allreduce(MPI_IN_PLACE, global_matrix, stop, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
		rmg_timings (DIAG_GLOB_SUMS, my_crtc () - time2);
		//QMD_dcopy (stop, global_matrix, ione, distAij, ione);
		cublasSetVector(num_states * num_states, sizeof( rmg_double_t ), global_matrix, ione, gpuAij, ione );

		// Now deal with the S operator
		time3 = my_crtc ();
		alpha = ct.vel;

		cublasSetVector(pbasis * num_states, sizeof( rmg_double_t ), pct.ns, ione, ct.gpu_temp, ione );
		cublasDgemm(ct.cublas_handle, cu_transT, cu_transN, num_states, num_states, pbasis,
				&alpha, ct.gpu_states, pbasis,
				ct.gpu_temp, pbasis,
				&beta,  ct.gpu_global_matrix, num_states );
		cublasGetVector(num_states * num_states, sizeof( rmg_double_t ), ct.gpu_global_matrix, ione, global_matrix, ione );


		time2 = my_crtc ();
		rmg_timings (DIAG_DGEMM, time2 - time3);

		// Reduce Sij and store it on the GPU
		//global_sums (global_matrix, &stop, pct.grid_comm);
		MPI_Allreduce(MPI_IN_PLACE, global_matrix, stop, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
		rmg_timings (DIAG_GLOB_SUMS, my_crtc () - time2);
		cublasSetVector(num_states * num_states, sizeof( rmg_double_t ), global_matrix, ione, gpuSij, ione );
		//QMD_dcopy (stop, global_matrix, ione, distSij, ione);

		/* Apply B operator on each wavefunction */
		time2 = my_crtc ();
//		subdiag_app_B (states, tmp_array2R);

		time3 = my_crtc ();
		alpha = ct.vel;
		rmg_timings (DIAG_APP_B, time3 - time2);

		cublasSetVector(pbasis * num_states, sizeof( rmg_double_t ), tmp_array2R, ione, ct.gpu_temp, ione );
		cublasDgemm(ct.cublas_handle, cu_transT, cu_transN, num_states, num_states, pbasis,
				&alpha, ct.gpu_states, pbasis,
				ct.gpu_temp, pbasis,
				&beta,  ct.gpu_global_matrix, num_states );
		cublasGetVector(num_states * num_states, sizeof( rmg_double_t ), ct.gpu_global_matrix, ione, global_matrix, ione );

		time2 = my_crtc ();
		rmg_timings (DIAG_DGEMM, time2 - time3);

		// Reduce Bij and leave in global_matrix
		//        global_sums (global_matrix, &stop, pct.grid_comm);
		MPI_Allreduce(MPI_IN_PLACE, global_matrix, stop, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
		rmg_timings (DIAG_GLOB_SUMS, my_crtc () - time2);



	}

	/*Create unitary matrix on GPU */
	for (st1 = 0; st1 < num_states*num_states; st1++)
	{
		distCij[st1] = 0.0;
	}
	for (st1 = 0; st1 < num_states;st1++) distCij[st1*num_states+st1] = 1.0;
	custat = cublasSetVector(num_states * num_states , sizeof( rmg_double_t ), distCij, ione, gpuCij, ione );

	time2 = my_crtc ();

	/*Get matrix that is inverse to B */
	{
		int *ipiv;

		my_calloc (ipiv, num_states, int);

		/* Transfer Bij which we left in global_matrix to the GPU. Inverse of B stored in gpuCij after dgesv call */
		cublasSetVector(num_states * num_states, sizeof( rmg_double_t ), global_matrix, ione, ct.gpu_global_matrix, ione );

		magma_dgesv_gpu( num_states, num_states, ct.gpu_global_matrix, num_states, ipiv, gpuCij, num_states, &info );
		// gpuCij holds B^-1

		if (info)
		{
			printf ("\n PE %d: p{d,z}gesv failed, info is %d", pct.gridpe, info);
			error_handler (" p{d,z}gesv failed");
		}

		my_free (ipiv);
	}

	/*Multiply inverse of B and and A */
	{
		char *trans = "n";

		/*B^-1*A */
		// ct.gpu_temp holds B^-1 and ct.gpu_work1 holds A
		cublasDgemm(ct.cublas_handle, cu_transN, cu_transN, num_states, num_states, num_states, &alpha1, 
				gpuCij, num_states, gpuAij, num_states, &beta1, gpuBij, 
				num_states );

		cublasDgemm(ct.cublas_handle, cu_transN, cu_transN, num_states, num_states, num_states, &alpha1, 
				gpuSij, num_states, gpuBij, num_states, &beta1, gpuCij, 
				num_states );

	}


	/****************** Find Matrix of Eigenvectors *****************************/
	/* Using lwork=-1, pdsyev should return minimum required size for the work array */
	{
		int eigs_found, eigvs_found, izero=0;
		int *iwork, *ifail, lwork;
		rmg_double_t *gap, lwork_tmp, *work2, qw1[10];
		int qw2[10];
		int liwork_tmp, liwork;

		my_malloc (ifail, num_states, int);
		lwork = -1;
		liwork = -1;

		// Get workspace reqs
		magma_dsyevd_gpu('V', 'L', num_states, gpuCij, num_states, eigs,
				distAij,  num_states,
				qw1, lwork,
				qw2, liwork,
				&info);

		lwork = (int)qw1[0];
		liwork = qw2[0];
		my_malloc (work2, lwork, rmg_double_t);
		my_malloc (iwork, liwork, int);
		//            info = rmg_dsygvd_gpu(num_states, gpuCij, num_states, gpuSij, num_states,
		//                   eigs, work2, lwork, iwork, liwork, distAij);
                if(ct.subdiag_driver == SUBDIAG_MAGMAFS) {
		    info = rmg_folded_spectrum_gpu(num_states, gpuCij, num_states, gpuSij, num_states,
				eigs, ct.gpu_host_work, lwork, iwork, liwork, distAij);
                }
                else {
		    info = rmg_dsygvd_gpu(num_states, gpuCij, num_states, gpuSij, num_states,
				eigs, ct.gpu_host_work, lwork, iwork, liwork, distAij);
                }

		if (info)
		{
			printf ("\n rmg_dsygvd_gpu failed, info is %d", info);
			error_handler ("rmg_dsygvd_gpu failed");
		}


		/*If subspace diagonalization is used everystep, use eigenvalues obtained here 
		 * as the correct eigenvalues*/
		if (ct.diag == 1)
			for (st1 = 0; st1 < ct.num_states; st1++) {
				states[st1].eig[0] = eigs[st1];
			}


		my_free (ifail);
		my_free (work2);
		my_free (iwork);

	}


	rmg_timings (DIAG_MATRIX_TIME, (my_crtc () - time2));

	/* Do the orbital update here. All data is already on the GPU */
	time2 = my_crtc ();

	alpha1 = 1.0;
	beta1 = 0.0;
	custat = cublasDgemm(ct.cublas_handle, cu_transN, cu_transN, pbasis, num_states, num_states,
			&alpha1, 
			ct.gpu_states, pbasis,
			gpuCij, num_states,
			&beta1,  ct.gpu_temp, pbasis );

	// Grab our data from the GPU.
	cublasGetVector(pbasis * num_states, sizeof( rmg_double_t ), ct.gpu_temp, ione, states->psiR, ione );


	rmg_timings (DIAG_WAVEUP_TIME, (my_crtc () - time2));


	/* release our temporary storage */
	my_free (work1R);
	my_free (eigs);
	my_free (vtot_eig);

	rmg_timings (DIAG_TIME, (my_crtc () - time1));

} // end subdiag_gamma_magma

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

// GPU specific version with itype=1,jobz=v,uplo=l
// assumes that a and b are already in gpu memory.
// Magma does not provide a routine that works as
// required so we put one together using magma routines
// and the cublas version of dtrsm.
int rmg_dsygvd_gpu(int n, rmg_double_t *a, int lda, rmg_double_t *b, int ldb, 
		rmg_double_t *w, rmg_double_t *work, int lwork, int *iwork, int liwork, rmg_double_t *wa)
{
	int ione=1, itype=1, info=0, idx;
	rmg_double_t rone = 1.0;
	cublasOperation_t cu_transT = CUBLAS_OP_T, cu_transN = CUBLAS_OP_N;
	cublasSideMode_t side=CUBLAS_SIDE_LEFT;
	cublasFillMode_t cuplo=CUBLAS_FILL_MODE_LOWER;
	cublasDiagType_t diag=CUBLAS_DIAG_NON_UNIT;

	//  Form a Cholesky factorization of B.
	//  This routine is buggy
	magma_dpotrf_gpu('L', n, b, ldb, &info);
	//cublasGetVector(n*n, sizeof( rmg_double_t ), b, ione, wa, ione );
	//dpotrf_("L", &n, wa, &n, &info);
	//cublasSetVector(n*n, sizeof( rmg_double_t ), wa, ione, b, ione );

	if( info != 0 ) {
		error_handler("dpotrf failure");
	}

	//  Transform problem to standard eigenvalue problem and solve.
	//   dsygst_( &itype, uplo, &n, a, &lda, b, &ldb, &info );
	//   dsyevd_( jobz, uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info );

	magma_dsygst_gpu(itype, 'L', n, a, lda, b, ldb, &info);
	if( info != 0 ) {
		error_handler("dsygst failure");
	}

	magma_dsyevd_gpu('V', 'L', n, a, lda, w,
			wa,  n,
			work, lwork,
			iwork, liwork,
			&info);

	if( info != 0 ) {
		error_handler("dsyevd failure");
	}

	//   For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
	//        backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
	//   dtrsm_( "Leftx", uplo, trans, "Non-unit", &n, &n, &rone, b, &ldb, a, &lda );
	//

	cublasDtrsm_v2 (ct.cublas_handle,side, cuplo, cu_transT, diag, n, n, &rone, b, ldb, a, lda);


	return 0;
} 
#else
// Empty stub 
void subdiag_gamma_magma (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc)
{
	error_handler("    This version of RMG was not built with Magma.\n");
}
#endif    // end #if over MAGMA_LIBS
#endif
#endif



#if GAMMA_PT
#if GPU_ENABLED
#if MAGMA_LIBS


#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


int rmg_folded_spectrum_gpu(int n, rmg_double_t *a, int lda, rmg_double_t *b, int ldb, 
		rmg_double_t *w, rmg_double_t *work, int lwork, int *iwork, int liwork, rmg_double_t *wa)
{

    int ione=1, itype=1, info=0, idx;
    rmg_double_t rone = 1.0;
    cublasOperation_t cu_transT = CUBLAS_OP_T, cu_transN = CUBLAS_OP_N;
    cublasSideMode_t side=CUBLAS_SIDE_LEFT;
    cublasFillMode_t cuplo=CUBLAS_FILL_MODE_LOWER;
    cublasDiagType_t diag=CUBLAS_DIAG_NON_UNIT;
    int ix, iy, eig_index, its, n_win, n_start, st, st1, i_width, ireps;
    char *trans="n";
    rmg_double_t alpha, beta=0.0, lambda, *d_p0, *d_p1, t1, t2;
    rmg_double_t *Vdiag, sum, *V, *G, *p0, *p1, *n_eigs, *tarr, *darr, *sarr, *Asave;
    rmg_double_t time1, time2, time3, r_width;
    int eig_step, eig_start, eig_stop, map, istride, ibase, omp_tid;
     
    // Folded spectrum method is parallelized over PE's. Each PE gets assigned
    // a subset of the eigenvectors.
    t1 = (rmg_double_t)n;
    t1 = t1 / ((rmg_double_t)NPES);
    t2 = t1 * (rmg_double_t)pct.gridpe;
    eig_start = (int)rint(t2);
    eig_stop = (int)rint(t1 + t2);
    eig_step = eig_stop - eig_start;
    if(pct.gridpe == (NPES - 1)) eig_stop = n;

    // Set width of window in terms of a percentage of n. Larger values will be slower but
    // exhibit behavior closer to full diagonalization.
    r_width = ct.folded_spectrum_width;
    t1 = (rmg_double_t)n;
    n_win = (int)(r_width * t1);

    // Find start of interval
    ix = n_win - eig_step;
    if(ix < 4)
        error_handler("Too few PE's to use folded spectrum method for this problem");
    if(ix % 2) {
        n_win++;
        ix = n_win - eig_step;
    }

    n_start = eig_start - ix/2;
    if(n_start < 0)n_start = 0;
    if((n_start + n_win) > n) {
        n_start = n - n_win;
    }

    my_malloc(Vdiag, n, rmg_double_t);
    my_malloc(p0, n, rmg_double_t);
    my_malloc(p1, n, rmg_double_t);
    n_eigs = distTij;

    if( cudaSuccess != cudaMalloc((void **)&d_p0 , n * sizeof(rmg_double_t) ))
        error_handler ("cudaMalloc failed for: d_p0\n");
    if( cudaSuccess != cudaMalloc((void **)&d_p1 , n * sizeof(rmg_double_t) ))
        error_handler ("cudaMalloc failed for: d_p1\n");

    time1=my_crtc();
    //  Form a Cholesky factorization of B.
    magma_dpotrf_gpu('L', n, b, ldb, &info);
    if( info != 0 ) {
        error_handler("dpotrf failure");
    }

    time2=my_crtc();
    Dprintf("DPOTRF1  = %12.6f",time2-time1);
    time1=my_crtc();

    //  Transform problem to standard eigenvalue problem
    magma_dsygst_gpu(itype, 'L', n, a, lda, b, ldb, &info);
    if( info != 0 ) {
        error_handler("dsygst failure");
    }

    time2=my_crtc();
    Dprintf("DSYGST  = %12.6f",time2-time1);
    time1=my_crtc();


    // We need to wait until a is diagonally dominant so we skip the first 3 steps
    if(ct.scf_steps < 3) {

        magma_dsyevd_gpu('V', 'L', n, a, lda, w,
                         wa,  n,
                         work, lwork,
                         iwork, liwork,
                         &info);
        if( info != 0 ) {
            error_handler("dsyevd failure");
        }
        //   For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
        //        backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
        //   dtrsm_( "Leftx", uplo, trans, "Non-unit", &n, &n, &rone, b, &ldb, a, &lda );
        //
        cublasDtrsm_v2 (ct.cublas_handle,side, cuplo, cu_transT, diag, n, n, &rone, b, ldb, a, lda);

    }
    else {
       
        V=tmp_arrayR;
        G=tmp_array2R;
        tarr = ct.gpu_host_temp1;
        Asave = ct.gpu_host_temp2;

        // AX=lambdaX  store a copy of A in Asave
        cublasGetVector(n * n, sizeof( rmg_double_t ), a, 1, Asave, 1 );

        // Zero out matrix of eigenvectors and eigenvalues
        for(idx=0;idx < n*n;idx++) V[idx] = 0.0;
        for(idx=0;idx < n;idx++) n_eigs[idx] = 0.0;
     

        time1=my_crtc();

        // Do the submatrix along the diagonal to get starting values for folded spectrum
        //--------------------------------------------------------------------
        for(ix = 0;ix < n_win;ix++){
            for(iy = 0;iy < n_win;iy++){
                G[ix*n_win + iy] = Asave[(n_start+ix)*n + n_start + iy];
            }
        }
        cublasSetVector(n_win * n_win, sizeof( rmg_double_t ), G, 1, a, 1 );
        magma_dsyevd_gpu('V', 'L', n_win, a, n_win, &w[n_start],
                        wa,  n_win,
                        work, lwork,
                        iwork, liwork,
                        &info);
        if( info != 0 ) {
                error_handler("dsyevd failure");
        }
        //--------------------------------------------------------------------

        cublasGetVector(n_win * n_win, sizeof( rmg_double_t ), a, 1, G, 1 );
        for(ix = 0;ix < n_win;ix++) {
            Vdiag[ix] = 1.0;
            if(G[ix*n_win + ix] < 0.0) Vdiag[ix] = -1.0;
        }

        // Store the eigen vector from the submatrix
        for(ix=0;ix<n_win;ix++) {
            if(((n_start+ix) >= eig_start) && ((n_start+ix) < eig_stop)) {
                for(iy=0;iy<n_win;iy++) {
                      V[(ix + n_start)*n + n_start + iy] = Vdiag[ix] * G[ix * n_win + iy];
                }
            }
        }


        time2=my_crtc();
        dprintf("SUBMATRIX = %12.6f",time2-time1);
        time1=my_crtc();

        // Apply folded spectrum to this PE's range of eigenvectors
        for(eig_index = eig_start;eig_index < eig_stop;eig_index++) {

                lambda = w[eig_index];
                n_eigs[eig_index] = lambda;

                for(ix=0;ix<ct.num_states;ix++){
                    Asave[ix*ct.num_states + ix] -= lambda;
                }
                alpha = 1.0;
                cublasSetVector( ct.num_states*ct.num_states, sizeof( rmg_double_t ), Asave, ione, ct.gpu_global_matrix, ione );

                // Restore matrix for next pass
                for(ix=0;ix<ct.num_states;ix++){
                    Asave[ix*ct.num_states + ix] += lambda;
                }

                cublasSetVector( n, sizeof( rmg_double_t ), &V[eig_index*n], ione, d_p0, ione );

                alpha = -0.5;
                beta = 0.0;
                for(its = 0;its < 8;its++) {
                    cublasDgemv_v2(ct.cublas_handle, cu_transN, n, n, &alpha, ct.gpu_global_matrix, n, d_p0, ione, &beta, d_p1, ione);
                    cublasDaxpy_v2(ct.cublas_handle, n, &rone, d_p1, ione, d_p0, ione);
                }
                // Renormalize
                cublasDnrm2(ct.cublas_handle, n, d_p0, ione, &t1);
                t1 = 1.0 / t1;
                cublasDscal(ct.cublas_handle, n, &t1, d_p0, ione);

            cublasGetVector(n, sizeof( rmg_double_t ), d_p0, ione, &V[eig_index*n], ione);
        }

        time2=my_crtc();
        dprintf("FOLDED SPECTRUM = %12.6f",time2-time1);
        time1=my_crtc();

        // Make sure all PE's have all eigenvectors. Possible optimization here would be to 
        // overlap computation in the above loop with communication here.
//        MPI_Allreduce(MPI_IN_PLACE, V, n*n, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        MPI_Allgatherv(MPI_IN_PLACE, eig_step*n, MPI_DOUBLE, V, fs_eigcounts, fs_eigstart, MPI_DOUBLE, pct.grid_comm);

        time2=my_crtc();
        dprintf("MPI_ALLREDUCE1  = %12.6f",time2-time1);
        time1=my_crtc();

        // Do the same for the eigenvalues
        // Could replace this with an MPI_Allgatherv as well
        MPI_Allreduce(MPI_IN_PLACE, n_eigs, n, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

        time2=my_crtc();
        dprintf("MPI_ALLREDUCE2  = %12.6f",time2-time1);

        // Copy summed eigs back to w
        for(idx = 0;idx < n;idx++) w[idx] = n_eigs[idx];

        time1=my_crtc();

        // Gram-Schmidt ortho for eigenvectors.
        alpha = 1.0;
        beta = 0.0;

        // Overlaps
        cublasSetVector( n * n, sizeof( rmg_double_t ), V, ione, a, ione );
        cublasDsyrk_v2 (ct.cublas_handle, cuplo, cu_transT, n, n, &alpha, a, n, &beta, ct.gpu_global_matrix, n);
        time2=my_crtc();
        dprintf("OVERLAPS  = %12.6f",time2-time1);
        time1=my_crtc();

        // Cholesky factorization
        magma_dpotrf_gpu('L', n, ct.gpu_global_matrix, n, &info);
        cublasGetVector( n * n, sizeof( rmg_double_t ), ct.gpu_global_matrix, ione, global_matrix, ione );

        time2=my_crtc();
        dprintf("CHOLESKY  = %12.6f",time2-time1);
        time1=my_crtc();

        // Get inverse of diagonal elements
        for(ix = 0;ix < n;ix++) tarr[ix] = 1.0 / global_matrix[n * ix + ix];

//----------------------------------------------------------------
        for(idx = 0;idx < n*n;idx++)G[idx] = 0.0;
#pragma omp parallel private(idx,st,st1,omp_tid,sarr)
{
        omp_tid = omp_get_thread_num();
        if(omp_tid == 0) my_malloc(darr, n * omp_get_num_threads(), rmg_double_t);
#pragma omp barrier

#pragma omp for schedule(static, 1) nowait
        for(idx = eig_start;idx < eig_stop;idx++) {

            sarr = &darr[omp_tid*n];

            for (st = 0; st < n; st++) sarr[st] = V[st*n + idx];

            for (st = 0; st < n; st++) {

                sarr[st] *= tarr[st];

                for (st1 = st+1; st1 < n; st1++) {
                    sarr[st1] -= global_matrix[st1 + n*st] * sarr[st];
                }

            }

            for (st = 0; st < n; st++) G[st*n + idx] = sarr[st];

        }
} // end omp section

        my_free(darr);

        time2=my_crtc();
        dprintf("GRAM  = %12.6f",time2-time1);
        time1=my_crtc();


        // A matrix transpose here would let us use an Allgatherv which would be
        // almost twice as fast for the communications part.
        MPI_Allreduce(MPI_IN_PLACE, G, n*n, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        cublasSetVector(n * n, sizeof( rmg_double_t ), G, 1, a, 1 );

        time2=my_crtc();
        dprintf("MPI_ALLREDUCE3  = %12.6f",time2-time1);
        time1=my_crtc();


        cublasDtrsm_v2 (ct.cublas_handle,side, cuplo, cu_transT, diag, n, n, &rone, b, ldb, a, lda);
        time2=my_crtc();
        dprintf("DTRSM  = %12.6f",time2-time1);

    }

    cudaFree(d_p1);
    cudaFree(d_p0);
    my_free(p1);
    my_free(p0);
    my_free(Vdiag);
    return 0;
} 
#endif    // end #if over MAGMA_LIBS
#endif
#endif

