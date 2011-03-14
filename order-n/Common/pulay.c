/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/* 
 *     Pulay  method subroutine 
 *
 *  	
 *
 * 	Input: 	step-- iteration step for Pulay,  MUST START from ZERO 
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



void pulay(int step, int N, double *xm, double *fm, int NsavedSteps, int preconditioning)
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
    double t1;
    double *x1, *x2, *f1, *f2;
    double *fi, *fj;
    int A_size;
    double scale = 0.5;
    double sd_step = 0.5;

    if (step == 0)
    {
        if (x == NULL)
            my_malloc( x, N * (NsavedSteps - 1), double );
        if (f == NULL)
            my_malloc( f, N * (NsavedSteps - 1), double );
        if (N > 0 && ((x == NULL) || (f == NULL)))
        {
            printf("pulay.c ---Could not allocate memory for x or f \n");
            fflush(NULL);
            exit(-1);
        }
    }
    if (step == 0)
    {
        scopy(&N, xm, &ione, x, &ione);
        scopy(&N, fm, &ione, f, &ione);

        if (preconditioning)
            precond(fm);
        alpha = -sd_step;
        saxpy(&N, &alpha, fm, &ione, xm, &ione);


    }
    else
    {
        size = step + 1;
        if (size > NsavedSteps)
            size = NsavedSteps;
        A_size = size + 1;
        s2 = A_size * A_size;

        for (i = 0; i < MAX_STEPS * MAX_STEPS; i++)
            A[i] = 0.0;
        for (i = 0; i < size; i++)
        {
            for (j = i; j < size; j++)
            {
                fi = f + i * N;
                fj = f + j * N;
                if (i == (size - 1))
                    fi = fm;
                if (j == (size - 1))
                    fj = fm;

                /*  get A(i,j)  */
                A[j * (size + 1) + i] = sdot(&N, fi, &ione, fj, &ione);
                A[i * (size + 1) + j] = A[j * (size + 1) + i];
            }

            b[i] = 0.0;
        }
        /*  Real_sum_all ( A, b )  if mutil-processing */
        global_sums(A, &s2);
        global_sums(b, &A_size);

        b[size] = 1.0;
        for (i = 0; i < size; i++)
        {
            A[i * (size + 1) + size] = 1.0;
            A[size * (size + 1) + i] = 1.0;
        }



        /*   b = A^(-1) * b     */
        sgesv(&A_size, &ione, A, &A_size, ipvt, b, &A_size, &info);

        if (pct.gridpe == 0)
        {
            printf("\n");
            for (i = 0; i < size; i++)
                printf("   pulay_b[%d]: %f  ", i, b[i]);
            printf("\n");
        }


        if (step <= (NsavedSteps - 2))
        {
            x1 = x + (size - 1) * N;
            f1 = f + (size - 1) * N;
            scopy(&N, xm, &ione, x1, &ione);
            scopy(&N, fm, &ione, f1, &ione);

            t1 = b[size - 1];
            sscal(&N, &t1, xm, &ione);
            for (i = 0; i < size - 1; i++)
            {
                x1 = x + i * N;
                saxpy(&N, &b[i], x1, &ione, xm, &ione);
            }

            t1 = 0.0;
            sscal(&N, &t1, fm, &ione);
            for (i = 0; i < size; i++)
            {
                t1 = -1.0 * b[i];
                f1 = f + i * N;
                saxpy(&N, &t1, f1, &ione, fm, &ione);
            }

            if (preconditioning)
                precond(fm);

            t1 = scale;
            saxpy(&N, &t1, fm, &ione, xm, &ione);

        }
        else
        {
            t1 = b[0];
            sscal(&N, &t1, x, &ione);
            for (i = 1; i < size - 1; i++)
            {
                x1 = x + i * N;
                saxpy(&N, &b[i], x1, &ione, x, &ione);
            }
            saxpy(&N, &b[size - 1], xm, &ione, x, &ione);

            t1 = -1.0 * b[0];
            sscal(&N, &t1, f, &ione);
            for (i = 1; i < size - 1; i++)
            {
                t1 = -1.0 * b[i];
                f1 = f + i * N;
                saxpy(&N, &t1, f1, &ione, f, &ione);
            }
            t1 = -1.0 * b[size - 1];
            saxpy(&N, &t1, fm, &ione, f, &ione);

            if (preconditioning)
                precond(f);

            t1 = scale;
            sscal(&N, &t1, f, &ione);
            saxpy(&N, &one, x, &ione, f, &ione);


            for (i = 0; i < size - 2; i++)
            {
                x1 = x + i * N;
                x2 = x + i * N + N;
                scopy(&N, x2, &ione, x1, &ione);
            }
            x1 = x + (size - 2) * N;
            scopy(&N, xm, &ione, x1, &ione);

            scopy(&N, f, &ione, xm, &ione);

            for (i = 0; i < size - 2; i++)
            {
                f1 = f + i * N;
                f2 = f + i * N + N;
                scopy(&N, f2, &ione, f1, &ione);
            }
            f1 = f + (size - 2) * N;
            scopy(&N, fm, &ione, f1, &ione);
        }

    }
}
