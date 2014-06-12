/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
                               get_dm_diag_p.c

    desca and ictxt are predefined in init_pe.c

	Requires the predefined parameters:
		NN:	size of the global matrix (my_scalapack.h)
		NPES:	number of processors (CFLAGS)

Documentation:

	LAPACK Working Note 94, A User's Guide to the BLACS v1.1

	Manpages: 

	WWW: http://www.netlib.org/scalapack/slug/scalapack_slug.html
	(ScaLapack user's guide)



			Fattebert J.-L., February 99 
			fatteber@nemo.physics.ncsu.edu                  

*/
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"



#include "my_scalapack.h"



void get_dm_diag_p(STATE * states, double *matS, double *X, double *hb)
{
	int num_states = ct.num_states;
    int ione = 1, izero = 0;    /* blas constants */
    char *uplo = "l", *jobz = "v";

    int info, numst = ct.num_states;
    double zero = 0., one = 1.;
    int st1, itype = 1;
    char diag, side, transa, transb;
    double *eigs;

    /* for parallel libraries */
    int nb, npcol, nprow;
    int mycol, myrow;
    int lwork;
    int  mxllda2;
    _fcd char_fcd1;
    _fcd char_fcd2;
    _fcd char_fcd3;


    void *RT = BeginRmgTimer("Diagonalization: scalapack");
    /* If I'm in the process grid, execute the program */
    if (pct.scalapack_myrow != -1)
    {


        /* 
         * SOLVE THE GENERALIZED EIGENVALUE PROBLEM:  m * z = lambda * matS * z 
         */

        /* Transform the generalized eigenvalue problem to a sStandard form */
        mxllda2 = MXLLDA * MXLCOL;
        scopy(&mxllda2, hb, &ione, uu_dis, &ione);
        scopy(&mxllda2, matS, &ione, l_s, &ione);

        char *range = "a";
        rmg_double_t vx = 0.0;
        rmg_double_t tol = 0.0;
        int eigs_found, eigvs_found;
        rmg_double_t orfac = -1.0;
        int *iwork, *ifail, *iclustr, lwork;
        rmg_double_t *gap, lwork_tmp, *work2;
        int liwork_tmp, liwork;

        my_malloc (ifail, num_states, int);
        my_malloc (eigs, num_states, double);
        my_malloc (iclustr, 2 * pct.scalapack_nprow * pct.scalapack_npcol, int);
        my_malloc (gap, pct.scalapack_nprow * pct.scalapack_npcol, rmg_double_t);

        lwork = -1;
        liwork = -1;


        PDSYGVX (&ione, jobz, range, uplo, &num_states, uu_dis, &ione, &ione, pct.desca,
                l_s, &ione, &ione, pct.desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
                &eigvs_found, eigs, &orfac, zz_dis, &ione, &ione, pct.desca, &lwork_tmp, &lwork,
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




        PDSYGVX (&ione, jobz, range, uplo, &num_states, uu_dis, &ione, &ione, pct.desca,
                l_s, &ione, &ione, pct.desca, &vx, &vx, &ione, &ione, &tol, &eigs_found,
                &eigvs_found, eigs, &orfac, zz_dis, &ione, &ione, pct.desca, work2, &lwork,
                iwork, &liwork, ifail, iclustr, gap, &info);


        if (info)
        {
            printf ("\n PDSYGVX failed, info is %d", info);
            error_handler ("PDSYGVX failed");
        }


        my_free (ifail);
        my_free (iclustr);
        my_free (gap);
        my_free (work2);
        my_free (iwork);

    }
    EndRmgTimer(RT);



    /* generalize the gamma for charge density matrix X */
    for (st1 = 0; st1 < numst; st1++)
    {
        work_matrix_row[st1] = states[st1].occupation[0];
    }

    diag_eig_matrix(gamma_dis, work_matrix_row, pct.desca);



    /* 
     * Build the density matrix X 
     */


    /* X = Z * gamma * Z^*  */
    side = 'r';
    uplo = 'l';
    char_fcd1 = &side;
    char_fcd2 = &uplo;
    PSSYMM(char_fcd1, char_fcd2, &numst, &numst, &one,
            gamma_dis, &ione, &ione, pct.desca,
            zz_dis, &ione, &ione, pct.desca, &zero, uu_dis, &ione, &ione, pct.desca);

    transa = 'n';
    transb = 't';
    char_fcd1 = &transa;
    char_fcd2 = &transb;
    PSGEMM(char_fcd1, char_fcd2, &numst, &numst, &numst, &one,
            uu_dis, &ione, &ione, pct.desca,
            zz_dis, &ione, &ione, pct.desca, &zero, X, &ione, &ione, pct.desca);


    my_barrier();

    MPI_Bcast(&eigs[0], numst, MPI_DOUBLE, 0, pct.grid_comm);
    for (st1 = 0; st1 < ct.num_states; st1++)
    {
        states[st1].eig[0] = eigs[st1];
    }

    /* my_free(work); */
    my_barrier();


}

