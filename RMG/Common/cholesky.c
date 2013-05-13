/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/cholesky.c *****
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
 *   void cholesky(rmg_double_t *a, int n)
 *   Computes the choleski decomposition of a symmetric positive 
 *   definite matrix if the matrix is not positive definite outputs
 *   an error message and halts program execution.
 * INPUT 
 *   a: The matrix whose decomposition is being computed.
 *   n: order of the matrix a(n,n)
 * OUTPUT
 *   a: the input matrix is overwritten
 * PARENTS
 *   subdiag_mpi.c subdiag_smp.c
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include "main.h"
#include <math.h>

#if 1
void cholesky(rmg_double_t * a, int n) {

    int info;
    char *uplo = "l";

    dpotrf(uplo, &n, a, &n, &info);

    if (info != 0)
        error_handler ("Matrix not positive definite or argument error.");

}

#else
void cholesky (rmg_double_t * a, int n)
{

    int i, j, k;
    rmg_double_t sum, *diagonal;

    my_malloc (diagonal, n, rmg_double_t);

    for (i = 0; i < n; i++)
    {
        for (j = i; j < n; j++)
        {

            for (sum = a[i * n + j], k = i - 1; k >= 0; k--)
                sum -= a[i * n + k] * a[j * n + k];
            if (i == j)
            {
                if (sum <= 0.0)
                    error_handler ("Matrix not positive definite");
                diagonal[i] = sqrt (sum);
            }
            else
            {
                a[j * n + i] = sum / diagonal[i];
                a[i * n + j] = a[j * n + i];
            }                   /* end if */

        }                       /* end for */
    }                           /* end for */


    for (i = 0; i < n; i++)
        a[i * n + i] = diagonal[i];


    my_free (diagonal);

}                               /* end cholesky */
#endif

/******/
