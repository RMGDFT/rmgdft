/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/* 
 *     KAIN  method subroutine 
 *
 *  	Krylov Subspace Accelerated Inexact Newton Method for  Linear and Nonlinear Equations 
 *
 * 	Input: 	step-- iteration step for KAIN,  MUST START from ZERO 
 * 		N  ---	Dinmention of array  xm and fm 	
 * 		xm ---	current sollution 
 * 		fm ---	current residual 
 * 		NsavedSteps ---	Number of steps saved, including the current step, Must greater than 2
 *
 * 	Output: xm ----- updated 
 * 
 *     	Note:  No D here, since preconditioning is done outside of this subroutine  
 *
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include "md.h"

#define 	MAX_STEPS 	300



void precond(double *x);



void kain(int step, int N, double *xm, double *fm, int NsavedSteps)
{
    static double *x;
    static double *f;
    double A[MAX_STEPS * MAX_STEPS];
    double b[MAX_STEPS];
    double c[MAX_STEPS];
    int ipvt[MAX_STEPS];
    int i, j;
    int ione = 1;
    int info;
    int size;
    int s2;
    char trans = 'N';
    double one = 1.0;
    double alpha;
    double beta;
    double gamma = -0.5;        /* KAIN  step length  */
    double sd_step = -0.1;
    double sum_ci = 0.0;
    double t1;
    double *x1, *x2, *f1, *f2;
    double xifj[MAX_STEPS * MAX_STEPS];
    double xmfj[MAX_STEPS];
    double xifm[MAX_STEPS];
    double xmfm;

    if (ct.scf_steps == 0)
    {
        my_malloc( x, N * (NsavedSteps - 1), double );
        my_malloc( f, N * (NsavedSteps - 1), double );
        if ((x == NULL) || (f == NULL))
        {
            printf("kain.c ---Could not allocate memory for x or f \n");
            fflush(NULL);
            exit(-1);
        }
    }
    if (step == 0)
    {
        scopy(&N, xm, &ione, x, &ione);
        scopy(&N, fm, &ione, f, &ione);

        precond(fm);

        saxpy(&N, &sd_step, fm, &ione, xm, &ione);
    }
    else
    {
        size = step;
        if (size > (NsavedSteps - 1))
            size = NsavedSteps - 1;

        for (i = 0; i < size; i++)
        {
            for (j = 0; j < size; j++)
            {
                x1 = x + i * N;
                f1 = f + j * N;
                xifj[i * size + j] = sdot(&N, x1, &ione, f1, &ione);
            }
            x1 = x + i * N;
            f1 = f + i * N;
            xmfj[i] = sdot(&N, xm, &ione, f1, &ione);
            xifm[i] = sdot(&N, x1, &ione, fm, &ione);
        }
        xmfm = sdot(&N, xm, &ione, fm, &ione);

        for (i = 0; i < size; i++)
        {
            for (j = 0; j < size; j++)
            {
                A[j * size + i] = xifj[i * size + j] - xmfj[j] - xifm[i] + xmfm;
            }
            b[i] = xmfm - xifm[i];
        }

        /*  Real_sum_all ( A, b )  if mutil-processing */
        s2 = size * size;
        global_sums(A, &s2);
        global_sums(b, &size);


        /*   b = A^(-1) * b     */
        sgesv(&size, &ione, A, &size, ipvt, b, &size, &info);


        if (pct.thispe == 0)
        {
            printf("\n");
            for (i = 0; i < size; i++)
            {
                sum_ci += b[i];
                printf(" KAIN_b[%d] = %f  ", i, b[i]);
            }
            printf("\n");
        }


        if (step <= (NsavedSteps - 2))
        {
            /*  xm = xm +sum(ci*(xi-xm)) - D*(sum(fj-fm)*cj+fm)   */
            x1 = x + size * N;
            f1 = f + size * N;
            scopy(&N, xm, &ione, x1, &ione);
            scopy(&N, fm, &ione, f1, &ione);

            t1 = 1.0 - sum_ci;
            sscal(&N, &t1, xm, &ione);
            for (i = 0; i < size; i++)
            {
                x1 = x + i * N;
                saxpy(&N, &b[i], x1, &ione, xm, &ione);
            }
            sscal(&N, &t1, fm, &ione);
            for (i = 0; i < size; i++)
            {
                f1 = f + i * N;
                saxpy(&N, &b[i], f1, &ione, fm, &ione);
            }

            precond(fm);
            saxpy(&N, &gamma, fm, &ione, xm, &ione);

        }
        else
        {
            t1 = b[0];
            sscal(&N, &t1, x, &ione);
            for (i = 1; i < size; i++)
            {
                x1 = x + i * N;
                saxpy(&N, &b[i], x1, &ione, x, &ione);
            }
            t1 = 1.0 - sum_ci;
            saxpy(&N, &t1, xm, &ione, x, &ione);

            t1 = b[0];
            sscal(&N, &t1, f, &ione);
            for (i = 1; i < size; i++)
            {
                f1 = f + i * N;
                saxpy(&N, &b[i], f1, &ione, f, &ione);
            }
            t1 = 1.0 - sum_ci;
            saxpy(&N, &t1, fm, &ione, f, &ione);

            precond(f);
            sscal(&N, &gamma, f, &ione);
            saxpy(&N, &one, x, &ione, f, &ione);

            for (i = 0; i < (size - 1); i++)
            {
                x1 = x + i * N;
                x2 = x + i * N + N;
                scopy(&N, x2, &ione, x1, &ione);
            }
            x1 = x + (size - 1) * N;
            scopy(&N, xm, &ione, x1, &ione);

            scopy(&N, f, &ione, xm, &ione);

            for (i = 0; i < (size - 1); i++)
            {
                f1 = f + i * N;
                f2 = f + i * N + N;
                scopy(&N, f2, &ione, f1, &ione);
            }
            f1 = f + (size - 1) * N;
            scopy(&N, fm, &ione, f1, &ione);

        }

    }
}
