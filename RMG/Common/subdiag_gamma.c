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
 *   void subdiag(STATE *states, REAL *vh, REAL *vnuc, REAL *vxc)
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
#include "main.h"

#if HYBRID_MODEL
#include "hybrid.h"
#include <pthread.h>
static pthread_attr_t diag_thread_attrs;
void subdiag_app_A_one_threaded(MG_THREAD_STRUCT *ss);
void subdiag_app_B_one_threaded(MG_THREAD_STRUCT *ss);
#endif

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
static void subdiag_app_A (STATE * states, REAL * a_psi, REAL * s_psi, REAL * vtot_eig);
static void subdiag_app_A_one (STATE * states, REAL * a_psi, REAL * s_psi, REAL * vtot_eig);
static void subdiag_app_B (STATE * states, REAL * b_psi);
static void subdiag_app_B_one (STATE * states, REAL * b_psi);
#else
void subdiag_app_A (STATE * states, REAL * a_psiR, REAL * a_psiI, REAL * s_psiR, REAL * s_psiI, REAL * vtot_eig);
static void subdiag_app_B (STATE * states, REAL * b_psiR, REAL * b_psiI);
#endif
#if GAMMA_PT
static void subdiag2_mpi (REAL * Aij, REAL * base_mem, REAL * tmp_psi);
#else
static void subdiag2_mpi (REAL * Aij, REAL * base_mem);
#endif
static void symmetrize_matrix (REAL * matrix, REAL * unity_matrix, int size, int *desca,
                               int local_size);
static void print_matrix2 (REAL * matrix, int size);
static void print_dist_matrix (REAL * dist_matrix, int global_size, int *desca);


void subdiag_gamma (STATE * states, REAL * vh, REAL * vnuc, REAL * vxc)
{
    int idx, st1;
    int num_states;
    int stop;
    int kidx;
    REAL *eigs;
    int ione = 1, izero = 0;    /* blas constants */
    char *uplo = "l", *jobz = "v";

    int info = 0;
    REAL time1, time2, time3;
    REAL *global_matrix;
    REAL *tmp_arrayR;
    REAL *tmp_array2R;
#if !GAMMA_PT
    REAL *tmp_arrayI;
#endif
    REAL *distAij, *distBij, *distCij, *distIij, *distSij;
    REAL *vtot, *vtot_eig;
	int dist_length, dist_stop, pbasis;

    pbasis = P0_BASIS;



    time1 = my_crtc ();

    if (pct.gridpe == 0)
        printf ("\n SUBSPACE DIAGONALIZATION");

    kidx = states[0].kidx;
    num_states = ct.num_states;


    /* Grab some temporary storage */
    stop = num_states * num_states;
#  if !GAMMA_PT
    stop = 2 * num_states * num_states;
#  endif

    /*Get memory for global matrices */

    my_calloc (global_matrix, stop, REAL);

    my_malloc (vtot_eig, P0_BASIS, REAL);
    my_malloc (eigs, num_states, REAL);


    /*Get vtot on coarse grid */
    my_malloc (vtot, FP0_BASIS, REAL);
    for (idx = 0; idx < FP0_BASIS; idx++)
        vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx];
    get_vtot_psi (vtot_eig, vtot, FG_NX);

    /*Release memory for vtot, do it already here since vtot is on fine grid */
    my_free (vtot);



    /*Temporary memory that will be used to calculate matrices and to update wavefunctions */
    my_malloc (tmp_arrayR, pbasis * ct.num_states, REAL);
    my_malloc (tmp_array2R, pbasis * ct.num_states, REAL);
#if !GAMMA_PT
    my_malloc (tmp_arrayI, pbasis * ct.num_states, REAL);
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


        /*Get memory for distributed matrices */
        my_calloc (distAij, dist_stop, REAL);
        my_calloc (distBij, dist_stop, REAL);
        my_calloc (distSij, dist_stop, REAL);
        my_calloc (distCij, dist_stop, REAL);
        my_calloc (distIij, dist_stop, REAL);
    }
    rmg_timings (DIAG_SCALAPACK_INIT, my_crtc () - time2);
    /********************* Scalapack should be initialized ******************************/








    /*Get matrix Aij */
    {
        char *trans = "t";
        char *trans2 = "n";
        REAL alpha = 1.0;
        REAL beta = 0.0;

        time2 = my_crtc ();

        /*Apply A operator on each wavefunction 
	 * S operator is also applied, the result is returned in tmp_array2R*/
        subdiag_app_A (states, tmp_arrayR, tmp_array2R, vtot_eig);

        time3 = my_crtc ();
        rmg_timings (DIAG_APP_A, time3 - time2);


        /*Global matrix will hold global A matrix */
        dgemm (trans, trans2, &num_states, &num_states, &pbasis, &alpha, states[0].psiR, &pbasis,
               tmp_arrayR, &pbasis, &beta, global_matrix, &num_states);


        time2 = my_crtc ();
        rmg_timings (DIAG_DGEMM, time2 - time3);


        /*Sum matrix over all processors */
        global_sums (global_matrix, &stop, pct.grid_comm);


        time3 = my_crtc ();
        rmg_timings (DIAG_GLOB_SUMS, time3 - time2);


        /*Distribute global matrix A */
        if (pct.scalapack_pe)
            distribute_mat (pct.desca, global_matrix, distAij, &num_states);

        time2 = my_crtc ();
        rmg_timings (DIAG_DISTMAT, time2 - time3);









        /* Now deal with S operator */



        time3 = my_crtc ();
        rmg_timings (DIAG_APP_S, time3 - time2);

        alpha = ct.vel;
        dgemm (trans, trans2, &num_states, &num_states, &pbasis, &alpha, states[0].psiR, &pbasis,
               tmp_array2R, &pbasis, &beta, global_matrix, &num_states);

        time2 = my_crtc ();
        rmg_timings (DIAG_DGEMM, time2 - time3);


        /*Sum matrix over all processors */
        global_sums (global_matrix, &stop, pct.grid_comm);

        time3 = my_crtc ();
        rmg_timings (DIAG_GLOB_SUMS, time3 - time2);


        /*Distribute global matrix A */
        if (pct.scalapack_pe)
            distribute_mat (pct.desca, global_matrix, distSij, &num_states);

        time2 = my_crtc ();
        rmg_timings (DIAG_DISTMAT, time2 - time3);









        /* Apply B operator on each wavefunction */
        subdiag_app_B (states, tmp_array2R);

        time3 = my_crtc ();
        rmg_timings (DIAG_APP_B, time3 - time2);

        dgemm (trans, trans2, &num_states, &num_states, &pbasis, &alpha, states[0].psiR, &pbasis,
               tmp_array2R, &pbasis, &beta, global_matrix, &num_states);


        time2 = my_crtc ();
        rmg_timings (DIAG_DGEMM, time2 - time3);



        /*Sum matrix over all processors */
        global_sums (global_matrix, &stop, pct.grid_comm);

        time3 = my_crtc ();
        rmg_timings (DIAG_GLOB_SUMS, time3 - time2);


        /*Distribute global matrix A */
        if (pct.scalapack_pe)
            distribute_mat (pct.desca, global_matrix, distBij, &num_states);

        time2 = my_crtc ();
        rmg_timings (DIAG_DISTMAT, time2 - time3);

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
        scopy (&dist_stop, distCij, &ione, distIij, &ione);


        /*Make Bij Hermitian */
        symmetrize_matrix (distBij, distIij, num_states, pct.desca, dist_length);




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
            REAL alpha[] = { 1.0, 0.0 };
            REAL beta[] = { 0.0, 0.0 };

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
            scopy (&dist_stop, distCij, &ione, distBij, &ione);
        }




        /*Make result symmetric (Hermitian) */
        symmetrize_matrix (distBij, distIij, num_states, pct.desca, dist_length);





#if 1
        /****************** Find Matrix of Eigenvectors *****************************/
        /* Using lwork=-1, pdsyev should return minimum required size for the work array */
        {
            char *range = "a";
            REAL vx = 0.0;
            REAL tol = 0.0;
            int eigs_found, eigvs_found;
            REAL orfac = 0.0;
            int *iwork, *ifail, *iclustr, lwork;
            REAL *gap, lwork_tmp, *work2;
            int liwork_tmp, liwork;

            my_malloc (ifail, num_states, int);
            my_malloc (iclustr, 2 * pct.scalapack_nprow * pct.scalapack_npcol, int);
            my_malloc (gap, pct.scalapack_nprow * pct.scalapack_npcol, REAL);
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

            my_malloc (work2, lwork, REAL);
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



        /*Free memory for distributed matrices */
        my_free (distSij);
        my_free (distIij);
        my_free (distCij);
        my_free (distBij);
        my_free (distAij);


    }                           /*end if (pct.scalapack_pe) */

    /*Other processors should have global_matrix set to 0 */
    else
        for (idx = 0; idx < stop; idx++)
            global_matrix[idx] = 0.0;


    rmg_timings (DIAG_MATRIX_TIME, (my_crtc () - time2));




    /*Finally, sum global_matrix over all PEs */
    time3 = my_crtc ();

    global_sums (global_matrix, &stop, pct.grid_comm);

    rmg_timings (DIAG_GLOB_SUMS, my_crtc () - time3);



    /*If some processors did not participate in Scalapack,
     * broadcast eigenvalues, since only Scalapack processors have updated eigenvalues*/
    if ((pct.scalapack_nprow * pct.scalapack_npcol != NPES) && (ct.diag == 1))
    {

        time2 = my_crtc ();

        MPI_Bcast (eigs, num_states, MPI_DOUBLE, 0, pct.grid_comm);

        /*Assign eigenvalues */
        for (st1 = 0; st1 < ct.num_states; st1++)
            states[st1].eig[0] = eigs[st1];

        rmg_timings (DIAG_BCAST_EIGS, (my_crtc () - time2));
    }




    time2 = my_crtc ();

    /* Do the orbital update in here */
#if GAMMA_PT
    subdiag2_mpi (global_matrix, states->psiR, tmp_arrayR);
#else
    subdiag2_mpi (global_matrix, states->psiR);
#endif


    rmg_timings (DIAG_WAVEUP_TIME, (my_crtc () - time2));







    /* release our temporary storage */

#if !GAMMA_PT
    my_free (tmp_arrayI);
#endif
    my_free (tmp_arrayR);
    my_free (tmp_array2R);
    my_free (global_matrix);
    my_free (vtot_eig);
    my_free (eigs);



    rmg_timings (DIAG_TIME, (my_crtc () - time1));

}


#  if GAMMA_PT
/*Applies A operator to all wavefunctions*/
static void subdiag_app_A (STATE * states, REAL * a_psi, REAL * s_psi, REAL * vtot_eig)
{
    int istate;
    STATE *sp;

#if HYBRID_MODEL
    {
        int st1, ist;
        pthread_t threads[THREADS_PER_NODE];
        MG_THREAD_STRUCT mst[THREADS_PER_NODE];
   
        pthread_attr_init( &diag_thread_attrs );
        pthread_attr_setschedpolicy( &diag_thread_attrs, SCHED_RR);
        scf_barrier_init(THREADS_PER_NODE);
        scf_tsd_init();
 
        // Each thread applies the operator to one wavefunction
        for(st1=0;st1 < ct.num_states;st1+=THREADS_PER_NODE) {
            for(ist = 0;ist < THREADS_PER_NODE;ist++) {
                mst[ist].sp = &states[st1 + ist];
                mst[ist].vtot = vtot_eig;
                mst[ist].tid = ist;
                mst[ist].p1 = &a_psi[(st1 + ist) * P0_BASIS];
                mst[ist].p2 = &s_psi[(st1 + ist) * P0_BASIS];
                pthread_create(&threads[ist], &diag_thread_attrs, (void *)subdiag_app_A_one_threaded, &mst[ist]);
            }

            for(ist = 0;ist < THREADS_PER_NODE;ist++) {
                pthread_join(threads[ist], NULL);
            }
        }

        scf_barrier_destroy();
        scf_tsd_delete();

    }
#else
    for (istate = 0; istate < ct.num_states; istate++) {
        sp = &states[istate];
        subdiag_app_A_one(sp, &a_psi[istate * P0_BASIS], &s_psi[istate * P0_BASIS], vtot_eig);
    }
#endif

}

#if HYBRID_MODEL
void subdiag_app_A_one_threaded(MG_THREAD_STRUCT *ss) {

    set_cpu_affinity();
    scf_tsd_set_value((void *)ss);
    scf_barrier_wait();
    subdiag_app_A_one(ss->sp, ss->p1, ss->p2, ss->vtot);

}

#endif


// Applies A operator to one wavefunction
static void subdiag_app_A_one (STATE *sp, REAL * a_psi, REAL * s_psi, REAL * vtot_eig)
{
    int kidx, idx, istate, sbasis;
    REAL *sg_twovpsi, *tmp_psi, *work2, *work1;
#    if MD_TIMERS
    REAL time1;
#    endif


    sbasis = sp->sbasis;
    my_malloc (work2, 2 * sbasis, REAL);
    sg_twovpsi = work2 + sbasis;
    kidx = 0;

    work1 = a_psi;

    tmp_psi = sp->psiR;


#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* Apply non-local operator to psi and store in work2 */
    app_nls (tmp_psi, NULL, work2, NULL, s_psi, NULL, ct.ions[0].newsintR, NULL, sp->istate, kidx);
#   if MD_TIMERS
    rmg_timings (DIAG_NL_TIME, (my_crtc () - time1));
#   endif



#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* Generate 2*V*psi and store it in a smoothing grid and store in sg_twovpsi */
    genvpsi (tmp_psi, sg_twovpsi, vtot_eig, work2, NULL, 0.0, PX0_GRID, PY0_GRID, PZ0_GRID);

#   if MD_TIMERS
        rmg_timings (DIAG_GENVPSI_TIME, (my_crtc () - time1));
#   endif

#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* B operating on 2*V*psi stored in work1 */
    app_cir_sixth (sg_twovpsi, work1, PX0_GRID, PY0_GRID, PZ0_GRID);
#   if MD_TIMERS
       rmg_timings (DIAG_APPCIR_TIME, (my_crtc () - time1));
#   endif

    /* Pack psi into smoothing array */
    //pack_ptos (sg_psi, tmp_psi, PX0_GRID, PY0_GRID, PZ0_GRID);


#   if MD_TIMERS
        time1 = my_crtc ();
#   endif

    /* A operating on psi stored in work2 */
    app_cil_sixth (tmp_psi, work2, PX0_GRID, PY0_GRID, PZ0_GRID, sp->hxgrid,
                   sp->hygrid, sp->hzgrid);

#   if MD_TIMERS
        rmg_timings (DIAG_APPCIL_TIME, (my_crtc () - time1));
#   endif

    for (idx = 0; idx < P0_BASIS; idx++)
        work1[idx] = 0.5 * ct.vel * (work1[idx] - work2[idx]);

    my_free (work2);

}                               /* subdiag_app_A_one */





/*Applies B operator to all wavefunctions*/
/*On input b_psi contains s-operator applied to wavefunction*/
static void subdiag_app_B (STATE * states, REAL * b_psi)
{
    int istate;
    STATE *sp;

#if HYBRID_MODEL
    {
        int st1, ist, istate;
        pthread_t threads[THREADS_PER_NODE];
        MG_THREAD_STRUCT mst[THREADS_PER_NODE];

        pthread_attr_init( &diag_thread_attrs );
//        pthread_attr_setschedpolicy( &diag_thread_attrs, SCHED_RR);
        scf_barrier_init(THREADS_PER_NODE);
        scf_tsd_init();

        // Each thread applies the operator to one wavefunction
        for(st1=0;st1 < ct.num_states;st1+=THREADS_PER_NODE) {
            for(ist = 0;ist < THREADS_PER_NODE;ist++) {
                mst[ist].sp = &states[st1 + ist];
                mst[ist].tid = ist;
                mst[ist].p1 = &b_psi[(st1 + ist) * P0_BASIS];
                pthread_create(&threads[ist], &diag_thread_attrs, (void *)subdiag_app_B_one_threaded, &mst[ist]);
            }

            for(ist = 0;ist < THREADS_PER_NODE;ist++) {
                pthread_join(threads[ist], NULL);
            }
        }

        scf_barrier_destroy();
        scf_tsd_delete();

    }
#else
    for (istate = 0; istate < ct.num_states; istate++) {
        sp = &states[istate];
        subdiag_app_B_one(sp, &b_psi[istate * P0_BASIS]);
    }
#endif

}


// Applies B operator to one wavefunction
static void subdiag_app_B_one (STATE *sp, REAL * b_psi)
{
    int istate, pbasis, ione=1;
    REAL *work2, *work1;
#    if MD_TIMERS
    REAL time1;
#    endif

    pbasis = sp->pbasis;

    my_malloc (work2, pbasis, REAL);

    work1 = b_psi;

    /*Pack S|psi> into smoothing array */
    //pack_ptos (sg_psi, work1, PX0_GRID, PY0_GRID, PZ0_GRID);
    scopy (&pbasis, work1, &ione, work2, &ione);


#   if MD_TIMERS
        time1 = my_crtc ();
#   endif
    /*B operating on S|psi> and store in work3 */
    app_cir_sixth (work2, work1, PX0_GRID, PY0_GRID, PZ0_GRID);
#   if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME2, (my_crtc () - time1));
#   endif

    my_free (work2);

}                               /* subdiag_app_B_one */


#if HYBRID_MODEL
void subdiag_app_B_one_threaded(MG_THREAD_STRUCT *ss) {

    set_cpu_affinity();
    scf_tsd_set_value((void *)ss);
    scf_barrier_wait();
    subdiag_app_B_one(ss->sp, ss->p1);

}

#endif



/* This routine is used to do the subspace rotation of the orbitals. Each
 * thread handles a specific portion of the real space domain.
 */
static void subdiag2_mpi (REAL * Aij, REAL * base_mem, REAL * tmp_psi)
{
    int idx;

    char *trans = "n";
    REAL alpha = 1.0;
    REAL beta = 0.0;
    int pbasis = P0_BASIS;
    int num_states = ct.num_states;


    dgemm (trans, trans, &pbasis, &num_states, &num_states, &alpha, base_mem, &pbasis, Aij,
           &num_states, &beta, tmp_psi, &pbasis);

    for (idx = 0; idx < num_states * pbasis; idx++)
        base_mem[idx] = tmp_psi[idx];


}                               /* end subdiag2_mpi */

#  else




/*Applies A operator to all wavefunctions*/
void subdiag_app_A (STATE * states, REAL * a_psiR, REAL * a_psiI, REAL * s_psiR, REAL * s_psiI, REAL * vtot_eig)
{
    int kidx, idx, istate, sbasis;
    REAL *sg_twovpsiR, *sg_twovpsiI, *tmp_psiR, *tmp_psiI, *work2R, *work2I,
        *work1R, *work1I;
    REAL *gx, *gy, *gz, *kdr;
    STATE *sp;
#    if MD_TIMERS
    REAL time1;
#    endif



    sbasis = states[0].sbasis;
    my_malloc (work2R, 8 * sbasis, REAL);
    work2I = work2R + sbasis;
    sg_twovpsiR = work2I + sbasis;
    sg_twovpsiI = sg_twovpsiR + sbasis;
    gx = sg_twovpsiI + sbasis;
    gy = gx + sbasis;
    gz = gy + sbasis;
    kdr = gz + sbasis;



    kidx = states[0].kidx;

    work1R = a_psiR;
    work1I = a_psiI;




    for (istate = 0; istate < ct.num_states; istate++)
    {

        sp = &states[istate];
        tmp_psiR = sp->psiR;
        tmp_psiI = sp->psiI;


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* Apply non-local operator to psi and store in work2 */
        app_nls (tmp_psiR, tmp_psiI, work2R, work2I, s_psiR, s_psiI, ct.ions[0].newsintR, ct.ions[0].newsintI, FALSE, kidx);
#    if MD_TIMERS
        rmg_timings (DIAG_NL_TIME, (my_crtc () - time1));
#    endif



        /* Pack psi into smoothing array and get image data */
        /*pack_ptos (sg_psiR, tmp_psiR, PX0_GRID, PY0_GRID, PZ0_GRID);
        pack_ptos (sg_psiI, tmp_psiI, PX0_GRID, PY0_GRID, PZ0_GRID);*/


        /* Apply the gradient operator to psi */
        app_grad (tmp_psiI, (P0_GRID *) gx, (P0_GRID *) gy, (P0_GRID *) gz);
        for (idx = 0; idx < P0_BASIS; idx++)
            kdr[idx] = (ct.kp[sp->kidx].kvec[0] * gx[idx] +
                        ct.kp[sp->kidx].kvec[1] * gy[idx] + ct.kp[sp->kidx].kvec[2] * gz[idx]);




#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* Generate 2*V*psi and store in sg_twovpsi */
        genvpsi (tmp_psiR, sg_twovpsiR, vtot_eig, work2R, kdr, ct.kp[sp->kidx].kmag, PX0_GRID,
                 PY0_GRID, PZ0_GRID);
#    if MD_TIMERS
        rmg_timings (DIAG_GENVPSI_TIME, (my_crtc () - time1));
#    endif


        /* Apply the gradient operator to psi */
        app_grad (tmp_psiR, (P0_GRID *) gx, (P0_GRID *) gy, (P0_GRID *) gz);
        for (idx = 0; idx < P0_BASIS; idx++)
            kdr[idx] = -(ct.kp[sp->kidx].kvec[0] * gx[idx] +
                         ct.kp[sp->kidx].kvec[1] * gy[idx] + ct.kp[sp->kidx].kvec[2] * gz[idx]);


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* Generate 2 * V * psiI and store it in a smoothing grid and store in sg_twovpsiI */
        genvpsi (tmp_psiI, sg_twovpsiI, vtot_eig, work2I, kdr, ct.kp[sp->kidx].kmag, PX0_GRID,
                 PY0_GRID, PZ0_GRID);
#    if MD_TIMERS
        rmg_timings (DIAG_GENVPSI_TIME, (my_crtc () - time1));
#    endif




#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* B operating on 2*V*psi stored in work1 */
        app_cir_sixth (sg_twovpsiR, work1R, PX0_GRID, PY0_GRID, PZ0_GRID);
        app_cir_sixth (sg_twovpsiI, work1I, PX0_GRID, PY0_GRID, PZ0_GRID);
#    if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME, (my_crtc () - time1));
#    endif



#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /* A operating on psi stored in work2 */
        app_cil_sixth (tmp_psiR, work2R, PX0_GRID, PY0_GRID, PZ0_GRID, sp->hxgrid,
                       sp->hygrid, sp->hzgrid);
        app_cil_sixth (tmp_psiI, work2I, PX0_GRID, PY0_GRID, PZ0_GRID, sp->hxgrid,
                       sp->hygrid, sp->hzgrid);

#    if MD_TIMERS
        rmg_timings (DIAG_APPCIL_TIME, (my_crtc () - time1));
#    endif

        for (idx = 0; idx < P0_BASIS; idx++)
        {
            work1R[idx] = 0.5 * ct.vel * (work1R[idx] - work2R[idx]);
            work1I[idx] = 0.5 * ct.vel * (work1I[idx] - work2I[idx]);
        }




        work1R += P0_BASIS;
        work1I += P0_BASIS;
	s_psiR += P0_BASIS;
	s_psiI += P0_BASIS;
    }

    my_free (work2R);

}                               /* subdiag_app_A */








/*Applies B operator to all wavefunctions*/
/*On input b_psi contains s-operator applied to wavefunction*/
void subdiag_app_B (STATE * states, REAL * b_psiR, REAL * b_psiI)
{
    int istate, ione=1;
    REAL *work2R, *work2I, *work1R, *work1I;
#    if MD_TIMERS
    REAL time1;
#    endif



    my_malloc (work2R, 2 * states[0].sbasis, REAL);
    work2I = work2R + states[0].sbasis;

    work1R = b_psiR;
    work1I = b_psiI;


    for (istate = 0; istate < ct.num_states; istate++)
    {

        /*Pack S|psi> into smoothing array */
        //pack_ptos (sg_psiR, work1R, PX0_GRID, PY0_GRID, PZ0_GRID);
        //pack_ptos (sg_psiI, work1I, PX0_GRID, PY0_GRID, PZ0_GRID);
	scopy (&states[0].sbasis, work1R, &ione, work2R, &ione);
	scopy (&states[0].sbasis, work1I, &ione, work2I, &ione);


#    if MD_TIMERS
        time1 = my_crtc ();
#    endif
        /*B operating on S|psi> and store in work3 */
        app_cir_sixth (work2R, work1R, PX0_GRID, PY0_GRID, PZ0_GRID);
        app_cir_sixth (work2I, work1I, PX0_GRID, PY0_GRID, PZ0_GRID);
#    if MD_TIMERS
        rmg_timings (DIAG_APPCIR_TIME2, (my_crtc () - time1));
#    endif


        work1R += P0_BASIS;
        work1I += P0_BASIS;
    }

    my_free (work2R);

}                               /* subdiag_app_B */








/* This routine is used to do the subspace rotation of the orbitals. Each
 * thread handles a specific portion of the real space domain.
 */
void subdiag2_mpi (REAL * Aij, REAL * base_mem)
{
    int idx, st1, st2;
    REAL *rptr;
    REAL *work1R, *work2R;
    REAL *work1I, *work2I;

    my_malloc (work1R, 4 * ct.num_states, REAL);
    work2R = work1R + ct.num_states;
    work1I = work2R + ct.num_states;
    work2I = work1I + ct.num_states;

    rptr = base_mem;

    for (idx = 0; idx < P0_BASIS; idx++)
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

static void print_matrix2 (REAL * matrix, int size)
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

static void print_dist_matrix (REAL * dist_matrix, int global_size, int *desca)
{
    REAL *glob_matrix;
    int stop;

    stop = global_size * global_size;

    my_calloc (glob_matrix, stop, REAL);


    if (pct.scalapack_pe)
        matgather (dist_matrix, desca, glob_matrix, global_size);


    /*Sum Aij over all PEs */
    global_sums (glob_matrix, &stop, pct.grid_comm);

    if (!pct.gridpe)
        print_matrix2 (glob_matrix, global_size);

    my_free (glob_matrix);


}



/*This works with distributed matrices*/
static void symmetrize_matrix (REAL * matrix, REAL * unity_matrix, int size, int *desca,
                               int local_size)
{
    int stop, ione = 1;
    REAL *temp_unity_matrix, *temp_matrix;
    REAL alpha[] = { 0.5, 0.0 };
    REAL beta[] = { 0.5, 0.0 };
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
    my_malloc (temp_matrix, stop, REAL);
    my_calloc (temp_unity_matrix, stop, REAL);

    /*Copy matrix into temp_matrix */
    scopy (&stop, matrix, &ione, temp_matrix, &ione);

    /*Local copy of unity matrix, this is done so that the unitary matrix that was passed here does not change */
    scopy (&stop, unity_matrix, &ione, temp_unity_matrix, &ione);


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
