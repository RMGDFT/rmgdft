/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 * 
 * INPUTS
 * 
 * OUTPUT
 *   nothing
 * PARENTS
 * 
 * CHILDREN
 * 
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "main.h"
#include "init_var.h"


void get_cholesky_real(double *matS)
{
    int ione = 1, numst = ct.num_states;
    int maxst;
    double time1, time2, time3;
    char uplo = 'l';
    int info, lwork;
    double rcond, anorm = 1.;
    double *work;
    int *iwork, liwork, locr;
    int nb, nproc, iproc;
    int mxllda, mxllda2;
    _fcd char_fcd1;

    int ix, iy, idx;

    nb = ct.scalapack_block_factor;
    maxst = ct.num_states;
    time3 = my_crtc();
#if DEBUG
    if (pct.gridpe == 0)
    {
        printf("\n Compute Cholesky decomposition \n");
    }
#endif


    mxllda = MXLLDA;
    mxllda2 = MXLLDA * MXLCOL;

    /* If I'm in the process grid, execute the program */
    if (pct.myrow != -1)
    {


        /*  Distribute matrices on the process grid */
        scopy(&mxllda2, matS, &ione, l_s, &ione);


        /* Compute the Cholesky decomposition of statearray */
        char_fcd1 = &uplo;
        PSPOTRF(char_fcd1, &numst, l_s, &ione, &ione, pct.desca, &info);
        if (info != 0)
        {
            printf(" PSPOTRF, output info=%d\n", info);
            fflush(NULL);
            exit(0);
        }

        /* Compute the conditioning of statearray */
        time1 = my_crtc();
        nproc = pct.nprow * pct.npcol;
        iproc = pct.gridpe;
        locr = ((numst / nb + 1) / nproc + 1) * nb + nb;
        lwork = locr * 5 + nb;
        liwork = locr;
        liwork *= 20;
        lwork *= 20;
        my_malloc( iwork, (size_t) (liwork), int );

        work = work_memory;
        PSPOCON(char_fcd1, &numst, l_s, &ione, &ione, pct.desca,
                &anorm, &rcond, work, &lwork, iwork, &liwork, &info);
        if (info != 0)
        {
		printf("\n work space %f %d \n", work[0], lwork);
            printf(" error in PSPOCON, info = %d\n", info);
            fflush(NULL);
            exit(0);
        }


        my_free(iwork);

/*
 *     RELEASE THE PROCESS GRID
 *     Free the BLACS context
 */
    }


    if (pct.gridpe == 0)
    {
        printf("\n Reciprocal of the condition number of the overlap matrix = %f\n", rcond);
        if (rcond < 0.00001)
        {
            printf(" Reciprocal of the condition number of the overlap matrix too small!!\n");
            fflush(NULL);
            exit(0);
        }
    }
    if (pct.gridpe == 0)
    {
        time2 = my_crtc();
        rmg_timings(COND_S_TIME, (time2 - time1));
        rmg_timings(CHOLESKY_TIME, (time2 - time3));
    }

}

/********/
