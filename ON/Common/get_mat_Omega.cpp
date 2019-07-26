/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

  This subroutine is used to calculate Omega = C * gamma * eigenvalue * C^T

*/
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "Scalapack.h"




void get_mat_Omega(STATE * states, double Omega[])
{
    int st1;
    int ione = 1, numst = ct.num_states;
    double zero = 0., one = 1.;
    char side, uplo = 'l', transa, transb;

    /* for parallel libraries */
    int myrow;
    char * char_fcd1;
    char * char_fcd2;


    myrow = pct.scalapack_myrow;

    /* If I'm in the process grid, execute the program */
    if (myrow != -1)
    {


        /* generilize the omega for mat_omega */
        for (st1 = 0; st1 < numst; st1++)
        {
            work_matrix_row[st1] = states[st1].occupation[0] * states[st1].eig[0];
        }

        diag_eig_matrix(gamma_dis, work_matrix_row, pct.desca);


        /* 
           Build the density matrix Omega 
           Omega = Z * gamma * eigenvalue * Z^*  
         */

        side = 'r';
        uplo = 'l';
        char_fcd1 = &side;
        char_fcd2 = &uplo;
        pdsymm (char_fcd1, char_fcd2, &numst, &numst, &one,
               gamma_dis, &ione, &ione, pct.desca,
               zz_dis, &ione, &ione, pct.desca, &zero, uu_dis, &ione, &ione, pct.desca);

        transa = 'n';
        transb = 't';
        char_fcd1 = &transa;
        char_fcd2 = &transb;
        pdgemm (char_fcd1, char_fcd2, &numst, &numst, &numst, &one,
               uu_dis, &ione, &ione, pct.desca,
               zz_dis, &ione, &ione, pct.desca, &zero, Omega, &ione, &ione, pct.desca);


    }

}

