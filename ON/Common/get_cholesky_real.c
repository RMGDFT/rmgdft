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
#include "prototypes_on.h"
#include "init_var.h"


void get_cholesky_real(double *matS)
{
    int ione = 1, numst = ct.num_states;
    char uplo = 'l';
    int info, lwork;
    double rcond, anorm = 1.;
    double *work;
    int *iwork, liwork, locr;
    int nb, nproc;
    int mxllda2;
    char *char_fcd1;


    nb = ct.scalapack_block_factor;
#if DEBUG
    if (pct.gridpe == 0)
    {
        printf("\n Compute Cholesky decomposition \n");
    }
#endif


    mxllda2 = MXLLDA * MXLCOL;

    /* If I'm in the process grid, execute the program */
    if (pct.scalapack_myrow != -1)
    {


        /*  Distribute matrices on the process grid */
        dcopy(&mxllda2, matS, &ione, l_s, &ione);


        /* Compute the Cholesky decomposition of statearray */
        char_fcd1 = &uplo;
        pdpotrf(char_fcd1, &numst, l_s, &ione, &ione, pct.desca, &info);
        if (info != 0)
        {
            printf(" pdpotrf, output info=%d\n", info);
            fflush(NULL);
            exit(0);
        }

        /* Compute the conditioning of statearray */
        nproc = pct.scalapack_nprow * pct.scalapack_npcol;
        locr = ((numst / nb + 1) / nproc + 1) * nb + nb;
        lwork = locr * 5 + nb;
        liwork = locr;
        liwork *= 20;
        lwork *= 20;
        my_malloc( iwork, (size_t) (liwork), int );

        work = work_memory;
        pdpocon(char_fcd1, &numst, l_s, &ione, &ione, pct.desca,
                &anorm, &rcond, work, &lwork, iwork, &liwork, &info);
        if (info != 0)
        {
	    printf("\n work space %f %d \n", work[0], lwork);
            printf(" error in pdpocon, info = %d\n", info);
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

}

/********/
