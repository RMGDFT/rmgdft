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

/* This subspace diagonalization function uses Scalapack libraries  */

/* This is used only for non-gamma point */
/* The reason we have two separate versions is that gamma point version has improvements
 * (matrix multiplications) when setting up matrices Aij, Bij, Cij and these improvements are not 
 * easily ported to non-gamma code, main problem is that matrix operations require that real 
 * and complex wavefunctions are continuous arrays, which they are not - we have 
 * states[0].psiR, states[0].psiI, states[1].psiR, states[1].psiI  etc
 * instead of required
 * states[0].psiR,  states[1].psiR, .... states[0].psiI, states[1].psiI
 * Chnaging this memory layout requires more time than I am willing to spend at the moment. 
 * Non-gamma point calculation is used only for small systems, and for small systems improvements to
 * subspace diagonalization  are not important since diagonalization takes only small percentage of whole
 * time.*/


#if !GAMMA_PT
static void subdiag1_mpi (int istate, STATE * states, REAL * Aij, REAL * Bij, REAL * vtot_eig);
static void subdiag2_mpi (REAL * Aij, REAL * base_mem);
static void symmetrize_matrix (REAL * matrix, REAL * unity_matrix, int size, int *desca,
                               int local_size);
static void print_matrix2 (REAL * matrix, int size);
static void print_dist_matrix (REAL * dist_matrix, int global_size, int *desca);


void subdiag_nongamma (STATE * states, REAL * vh, REAL * vnuc, REAL * vxc)
{
    int idx, st1;
	int num_states;
    int stop;
    int kidx;
    REAL *work1;
    int ione = 1, izero = 0;    /* blas constants */
    char *uplo = "l", *jobz = "v";

    int info = 0;
    REAL time1, time2, time3;
    REAL *Aij, *Bij, *Cij;
    REAL *distAij, *distBij, *distCij, *distIij;
    REAL *vtot, *vtot_eig;
    int dist_length, dist_stop;


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
    my_malloc (Aij, stop, REAL);
    my_malloc (Bij, stop, REAL);
    my_calloc (Cij, stop, REAL);

    my_malloc (work1, 8 * num_states, REAL);
    my_malloc (vtot_eig, P0_BASIS, REAL);


    /*Get vtot on coarse grid */
    my_malloc (vtot, FP0_BASIS, REAL);
    for (idx = 0; idx < FP0_BASIS; idx++)
        vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx];
    get_vtot_psi (vtot_eig, vtot, FG_NX);

    /*Release memory for vtot, do it already here since vtot is on fine grid */
    my_free (vtot);





    time2 = my_crtc ();
    /*This setups matrices Aij and Bij
     * Aij is <psi|H|psi>, Bij is <psi|B|psi>*/
    for (st1 = 0; st1 < num_states; st1++)
        subdiag1_mpi (st1, states, Aij, Bij, vtot_eig);

    rmg_timings (DIAG_SUBDIAG1_TIME, (my_crtc () - time2), 0);


    /* Sum A and B overlap matrices over all processors */
    time3 = my_crtc ();

    global_sums (Bij, &stop, pct.grid_comm);
    global_sums (Aij, &stop, pct.grid_comm);

    rmg_timings (DIAG_GLOB_SUMS, my_crtc () - time3, 0);



    /*Setup Cij to be unitary matrix */
    for (st1 = 0; st1 < num_states; st1++)
    {
#if GAMMA_PT
        Cij[st1 * num_states + st1] = 1.0;
#else
        Cij[2 * (st1 * num_states + st1)] = 1.0;
#endif
    }



#if 0
    if (pct.gridpe == 0)
    {
        printf ("\n\n Bij");
        print_matrix2 (Bij, num_states);

        printf ("\n\n Aij");
        print_matrix2 (Aij, num_states);

    }
#endif


    time2 = my_crtc ();









    if (pct.scalapack_pe)
    {

        /*Length of distributed matrices (different on each processor) */
        dist_length =
            NUMROC (&num_states, &pct.desca[4], &pct.scalapack_myrow, &izero,
                    &pct.scalapack_nprow) * NUMROC (&num_states, &pct.desca[4], &pct.scalapack_mycol,
                                                    &izero, &pct.scalapack_npcol);

        /* Make sure that we do not ask for memory of size 0 or thing like that */
        if (dist_length == 0)
            dist_length = 1;

        /*This holds number of doubles on each PE */
        dist_stop = dist_length;
#if !GAMMA_PT
        dist_stop *= 2;
#endif




        /*Get memory for distributed matrices */
        my_calloc (distAij, dist_stop, REAL);
        my_calloc (distBij, dist_stop, REAL);
        my_calloc (distCij, dist_stop, REAL);
        my_calloc (distIij, dist_stop, REAL);




        /*Distribute matrices */
        distribute_mat (pct.desca, Aij, distAij, &num_states);
        distribute_mat (pct.desca, Bij, distBij, &num_states);
        distribute_mat (pct.desca, Cij, distCij, &num_states);


        /*keep an extra copy of distributed unitary matrix */
        scopy (&dist_stop, distCij, &ione, distIij, &ione);


        /*Make Bij Hermitian */
#if !GAMMA_PT
        symmetrize_matrix (distBij, distIij, num_states, pct.desca, dist_length);
#endif



#if 0
        if (pct.gridpe == 0)
            printf ("\n\n Distributed Bij");
        print_dist_matrix (distBij, num_states, pct.desca);

        if (pct.gridpe == 0)
            printf ("\n\n Distributed Aij");
        print_dist_matrix (distAij, num_states, pct.desca);
#endif


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
        }




        /*Make result symmetric (Hermitian) */
        symmetrize_matrix (distBij, distIij, num_states, pct.desca, dist_length);





        /****************** Find Matrix of Eigenvectors *****************************/
        /* Using lwork=-1, pdsyev should return minimum required size for the work array */
        {
            int lwork = -1, lrwork = -1;
            REAL t1[2], t2[2], *work2, *work3;



#if GAMMA_PT
            PDSYEV (jobz, uplo, &num_states, distBij, &ione, &ione, pct.desca, work1, distAij, &ione,
                    &ione, pct.desca, t1, &lwork, &info);
#else
            PCHEEV (jobz, uplo, &num_states, distBij, &ione, &ione, pct.desca, work1, distAij, &ione,
                    &ione, pct.desca, t1, &lwork, t2, &lrwork, &info);
            lrwork = t2[0] + 1;
#endif

            /*Adding 1  to take care of case when t1 is, say 14.99, then lwork would only be 14 */
            lwork = t1[0] + 1;



#if GAMMA_PT
            my_malloc (work2, lwork, REAL);
            PDSYEV (jobz, uplo, &num_states, distBij, &ione, &ione, pct.desca, work1, distAij, &ione,
                    &ione, pct.desca, work2, &lwork, &info);
#else
            /*Calculate lwork first */
            my_malloc (work2, 2 * lwork, REAL);
            my_malloc (work3, 2 * lrwork, REAL);
            PCHEEV (jobz, uplo, &num_states, distBij, &ione, &ione, pct.desca, work1, distAij, &ione,
                    &ione, pct.desca, work2, &lwork, work3, &lrwork, &info);
#endif


            if (info)
                error_handler ("PDSYEV failed");

            my_free (work2);
#if !GAMMA_PT
            my_free (work3);
#endif
        }



        /*Gather result onto Aij */
        matgather (distAij, pct.desca, Aij, num_states);



        /*Free memory for distributed matrices */
        my_free (distIij);
        my_free (distCij);
        my_free (distBij);
        my_free (distAij);


    }                           /*end if (pct.scalapack_pe) */

    /*Other processors should have Aij set to 0 */
    else
        for (idx = 0; idx < stop; idx++)
            Aij[idx] = 0.0;



    /*Finally, sum Aij over all PEs */
    time3 = my_crtc ();

    global_sums (Aij, &stop, pct.grid_comm);

    rmg_timings (DIAG_GLOB_SUMS, my_crtc () - time3, 0);


#if 0
    if (pct.gridpe == 0)
    {
        printf ("\n\n Final final matrix");
        print_matrix2 (Aij, ct.num_states);

    }
#endif


    rmg_timings (DIAG_MATRIX_TIME, (my_crtc () - time2), 0);
    time2 = my_crtc ();


    /* Do the orbital update in here */
    subdiag2_mpi (Aij, states->psiR);


    rmg_timings (DIAG_WAVEUP_TIME, (my_crtc () - time2), 0);





    /* release our temporary storage */
    my_free (work1);
    my_free (Cij);
    my_free (Bij);
    my_free (Aij);
    my_free (vtot_eig);


    rmg_timings (DIAG_TIME, (my_crtc () - time1), 0);


    /*fflush(NULL);
       exit(0); */

}



static void subdiag1_mpi (int istate, STATE * states, REAL * Aij, REAL * Bij, REAL * vtot_eig)
{
    int idx, st2, ione = 1;
    int dimx, dimy, dimz, kidx;
    REAL *sg_twovpsiR, *sg_twovpsiI, *tmp_psiR, *tmp_psiI;
    REAL *work1R, *work1I;
    REAL *work2R, *work2I, *work3R, *work3I, *work4R, *work4I;
    REAL *gx, *gy, *gz, *kdr;
    REAL s1, s2, s3, s4;
    STATE *sp, *sp1;
#    if MD_TIMERS
    REAL time1;
#    endif

    sp = &states[istate];
    kidx = sp->kidx;
    my_malloc (work4R, 15 * sp->sbasis, REAL);
    work4I = work4R + sp->sbasis;
    work1R = work4I+ sp->sbasis;
    work1I = work1R + sp->sbasis;
    work2R = work1I + sp->sbasis;
    work2I = work2R + sp->sbasis;
    work3R = work2I + sp->sbasis;
    work3I = work3R + sp->sbasis;
    sg_twovpsiR = work3I + sp->sbasis;
    sg_twovpsiI = sg_twovpsiR + sp->sbasis;
    gx = sg_twovpsiI + sp->sbasis;
    gy = gx + sp->sbasis;
    gz = gy + sp->sbasis;
    kdr = gz + sp->sbasis;

    dimx = sp->dimx;
    dimy = sp->dimy;
    dimz = sp->dimz;

    /*No need to translate data from wavefunctions into tmp_psi, getting pointers should be enough */
    /*gather_psi(tmp_psiR, tmp_psiI, sp, 0); */
    tmp_psiR = sp->psiR;
    tmp_psiI = sp->psiI;


#    if MD_TIMERS
    time1 = my_crtc ();
#    endif

    /* Apply non-local operator to psi and store in work2 */
    app_nls (tmp_psiR, tmp_psiI, work2R, work2I, work3R, work3I, ct.ions[0].newsintR, ct.ions[0].newsintI, sp->istate, sp->kidx);


#    if MD_TIMERS
    rmg_timings (DIAG_NLS_TIME, (my_crtc () - time1), 0);
#    endif

    /*Pack work3 into smoothing array */
    /*pack_ptos (sg_psiR, work3R, dimx, dimy, dimz);
    pack_ptos (sg_psiI, work3I, dimx, dimy, dimz);*/

#    if MD_TIMERS
    time1 = my_crtc ();
#    endif
    /*B operating on S|psi> and store in work3 */
    app_cir_sixth (work3R, work4R, dimx, dimy, dimz);
    app_cir_sixth (work3I, work4I, dimx, dimy, dimz);

#    if MD_TIMERS
    rmg_timings (DIAG_APPCIR_TIME, (my_crtc () - time1), 0);
#    endif


    /* Pack psi into smoothing array and get image data */
    /*pack_ptos (sg_psiR, tmp_psiR, dimx, dimy, dimz);
    pack_ptos (sg_psiI, tmp_psiI, dimx, dimy, dimz);*/


    /* Apply the gradient operator to psi */
    app_grad (tmp_psiI, (P0_GRID *) gx, (P0_GRID *) gy, (P0_GRID *) gz);
    for (idx = 0; idx < P0_BASIS; idx++)
        kdr[idx] = (ct.kp[sp->kidx].kvec[0] * gx[idx] +
                    ct.kp[sp->kidx].kvec[1] * gy[idx] + ct.kp[sp->kidx].kvec[2] * gz[idx]);


    /* Generate 2 * V * psiR and store it in sg_twovpsiR */
    genvpsi (tmp_psiR, sg_twovpsiR, vtot_eig, work2R, kdr, ct.kp[sp->kidx].kmag, dimx, dimy, dimz);

    /* Apply the gradient operator to psi */
    app_grad (tmp_psiR, (P0_GRID *) gx, (P0_GRID *) gy, (P0_GRID *) gz);
    for (idx = 0; idx < P0_BASIS; idx++)
        kdr[idx] = -(ct.kp[sp->kidx].kvec[0] * gx[idx] +
                     ct.kp[sp->kidx].kvec[1] * gy[idx] + ct.kp[sp->kidx].kvec[2] * gz[idx]);


    /* Generate 2 * V * psiI and store it in a smoothing grid and store in sg_twovpsiI */
    genvpsi (tmp_psiI, sg_twovpsiI, vtot_eig, work2I, kdr, ct.kp[sp->kidx].kmag, dimx, dimy, dimz);

#    if MD_TIMERS
    time1 = my_crtc ();
#    endif

    /* B operating on 2*V*psi stored in work1 */
    app_cir_sixth (sg_twovpsiR, work1R, dimx, dimy, dimz);
    app_cir_sixth (sg_twovpsiI, work1I, dimx, dimy, dimz);

#    if MD_TIMERS
    rmg_timings (DIAG_APPCIR_TIME, (my_crtc () - time1), 0);
#    endif

#    if MD_TIMERS
    time1 = my_crtc ();
#    endif

    /* A operating on psi stored in work2 */
    app_cil_sixth (tmp_psiR, work2R, dimx, dimy, dimz,
                   sp->hxgrid, sp->hygrid, sp->hzgrid);
    app_cil_sixth (tmp_psiI, work2I, dimx, dimy, dimz,
                   sp->hxgrid, sp->hygrid, sp->hzgrid);

#    if MD_TIMERS
    rmg_timings (DIAG_APPCIL_TIME, (my_crtc () - time1), 0);
#    endif

    for (idx = 0; idx < sp->pbasis; idx++)
        work1R[idx] = 0.5 * ct.vel * (work1R[idx] - work2R[idx]);

    for (idx = 0; idx < sp->pbasis; idx++)
        work1I[idx] = 0.5 * ct.vel * (work1I[idx] - work2I[idx]);


    /* Compute the complex overlap matrices here */

#    if MD_TIMERS
    time1 = my_crtc ();
#    endif
    for (st2 = 0; st2 < ct.num_states; st2++)
    {

        sp1 = &states[st2];
        /*No need to translate data from wavefunctions into tmp_psi, getting pointers should be enough */
        /*gather_psi(tmp_psiR, tmp_psiI, sp1, 0); */
        tmp_psiR = sp1->psiR;
        tmp_psiI = sp1->psiI;

        s1 = QMD_sdot (sp->pbasis, tmp_psiR, ione, work1R, ione);
        s2 = QMD_sdot (sp->pbasis, tmp_psiR, ione, work1I, ione);
        s3 = QMD_sdot (sp->pbasis, tmp_psiI, ione, work1R, ione);
        s4 = QMD_sdot (sp->pbasis, tmp_psiI, ione, work1I, ione);

        Aij[2 * (sp->istate * ct.num_states + st2)] = (s1 + s4);
        Aij[2 * (sp->istate * ct.num_states + st2) + 1] = (s2 - s3);


        s1 = QMD_sdot (sp->pbasis, tmp_psiR, ione, work4R, ione);
        s2 = QMD_sdot (sp->pbasis, tmp_psiR, ione, work4I, ione);
        s3 = QMD_sdot (sp->pbasis, tmp_psiI, ione, work4R, ione);
        s4 = QMD_sdot (sp->pbasis, tmp_psiI, ione, work4I, ione);

        Bij[2 * (sp->istate * ct.num_states + st2)] = ct.vel * (s1 + s4);
        Bij[2 * (sp->istate * ct.num_states + st2) + 1] = ct.vel * (s2 - s3);


    }                           /* st2 */

    /*Bij should be hermitian matrix, this will be enforced later */

#    if MD_TIMERS
    rmg_timings (DIAG_SUBDIAG1_LOOP_TIME, (my_crtc () - time1), 0);
#    endif

    my_free (work4R);


}                               /* end subdiag1_mpi */



/* This routine is used to do the subspace rotation of the orbitals. Each
 * thread handles a specific portion of the real space domain.
 */
static void subdiag2_mpi (REAL * Aij, REAL * base_mem)
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



static void print_matrix2 (REAL * matrix, int size)
{
    int i, j;

    //printf("\n\n");

    for (i = 0; i < size; i++)
    {
        printf ("\n");
        for (j = 0; j < size; j++)
            printf ("%g  ", matrix[i * size + j]);
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
