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
#include "md.h"



#include "my_scalapack.h"



void get_dm_diag_p(STATE * states, double l_s[], double X[], double hb[])
{
    int i, ione = 1, info, numst = ct.num_states;
    double zero = 0., one = 1.;
    int st1, itype = 1;
    char diag, side, uplo = 'l', jobz = 'v', transa, transb, char_tmp;
    double time1, time2, time3;
    double val[NN];

    /* for parallel libraries */
    int nb = NB, npcol, nprow;
    int mycol, myrow;
    int lwork;
    int rsrc = 0, csrc = 0;
    int mxllda, mxllda2;
    double scale;
    double *work;
    _fcd char_fcd1;
    _fcd char_fcd2;
    _fcd char_fcd3;


    int icrow, iccol, mpc0, nqc0, sizemqrleft, qrmem, i1, itwo = 2, izero = 0;
    int npes = NPES, nn = NN, ldc, nrc;

    my_barrier();
    time3 = my_crtc();


    mxllda = MXLLDA;
    mxllda2 = MXLLDA * MXLLDA;
    npcol = pct.npcol;
    nprow = pct.nprow;

    mycol = pct.mycol;
    myrow = pct.myrow;

    icrow = INDXG2P(&itwo, &nb, &myrow, &izero, &nprow);
    iccol = INDXG2P(&ione, &nb, &mycol, &izero, &npcol);
    mpc0 = NUMROC(&nn, &nb, &myrow, &icrow, &nprow);
    nqc0 = NUMROC(&nn, &nb, &mycol, &iccol, &npcol);
    nrc = NUMROC(&nn, &nb, &myrow, &izero, &npes);
    ldc = max(1, nrc);
    sizemqrleft = max((NB * (NB - 1)) / 2, (nqc0 + mpc0) * NB) + NB * NB;
    sizemqrleft *= 2;
    qrmem = 2 * NN - 2;
    lwork = 5 * NN + NN * ldc + max(sizemqrleft, qrmem) + 1;


    /* my_malloc_init( work, lwork, REAL ); */
    work = work_memory;


    /* If I'm in the process grid, execute the program */
    if (myrow != -1)
    {


        /* 
         * SOLVE THE GENERALIZED EIGENVALUE PROBLEM:  m * z = lambda * matS * z 
         */
        time1 = my_crtc();

        /* Transform the generalized eigenvalue problem to a standard form */
        scopy(&mxllda2, hb, &ione, uu_dis, &ione);

        char_fcd1 = &uplo;

        PSSYGST(&itype, char_fcd1, &numst, uu_dis, &ione, &ione, pct.desca,
                l_s, &ione, &ione, pct.desca, &scale, &info);
        if (info != 0)
        {
            printf(" get_dm_diag_p: PSSYGST, output info=%d\n", info);
            fflush(NULL);
            exit(0);
        }
        if (myrow == 0 && mycol == 0)
        {
            time2 = my_crtc();
            rmg_timings(PSSYGST_TIME, (time2 - time1), 0);
        }


        /* solve a standard symmetric eigenvalue problem */
        time1 = my_crtc();

        char_fcd1 = &jobz;
        char_fcd3 = &uplo;


        PSSYEV(char_fcd1, char_fcd3, &numst, uu_dis, &ione, &ione, pct.desca,
               val, zz_dis, &ione, &ione, pct.desca, work, &lwork, &info);
        if (info != 0)
        {
            printf(" get_dm_diag_p: PSSYEV, output info=%d\n", info);
            fflush(NULL);
            exit(0);
        }
        if (myrow == 0 && mycol == 0)
        {
            time2 = my_crtc();
            rmg_timings(PSSYEV_TIME, (time2 - time1), 0);
        }


        /* generalize the gamma for charge density matrix X */
        for (st1 = 0; st1 < numst; st1++)
        {
            work_matrix_row[st1] = 0.5 * states[st1].occupation;
        }

        diag_eig_matrix(gamma_dis, work_matrix_row, pct.desca);


        /* Get the eigenvectors Z of the generalized eigenvalue problem */

        /* Solve Z=L**(-T)*U */
        time1 = my_crtc();
        diag = 'n';
        transa = 't';
        char_fcd1 = &uplo;
        char_fcd2 = &transa;
        char_fcd3 = &diag;
        PSTRTRS(char_fcd1, char_fcd2, char_fcd3, &numst, &numst,
                l_s, &ione, &ione, pct.desca, zz_dis, &ione, &ione, pct.desca, &info);
        if (info != 0)
        {
            printf(" get_dm_diag_p: PSTRTRS, output info=%d\n", info);
            fflush(NULL);
            exit(0);
        }
        if (myrow == 0 && mycol == 0)
        {
            time2 = my_crtc();
            rmg_timings(PSTRTRS_TIME, (time2 - time1), 0);
        }



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




    }


    my_barrier();

    MPI_Bcast(&val[0], numst, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (st1 = 0; st1 < ct.num_states; st1++)
    {
        states[st1].eig = val[st1];
    }

    /* my_free(work); */
    my_barrier();
    if (pct.gridpe == 0)
    {
        time2 = my_crtc();
        rmg_timings(DIAG_TIME, time2 - time3, 0);
    }


}

