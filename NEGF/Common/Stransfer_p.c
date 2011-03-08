/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
 *   Iterative construction of the transfer matrix, 
 *	as Lopez-Sancho & Rubio, J.Phys.F: Met. Phs., v.14, 1205 (1984) 
 *		and ibid. v.15, 851 (1985)
 *
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "md.h"
#include "pmo.h"

#define 	MAX_STEP 	100

static void my_print_sum (doublecomplex *matrix, int nmax, int index);

void Stransfer_p (doublecomplex * tot, doublecomplex * tott, REAL * H00, REAL * H01,
        REAL * S00, REAL * S01, REAL eneR, REAL eneI, int iprobe)
{

    REAL converge1, converge2;
    doublecomplex *tau, *taut, *tsum, *tsumt, *t11, *t12, *s1, *s2;
    doublecomplex *ch0, *ch1;
    doublecomplex alpha, beta;
    int info;
    int *ipiv;
    int i, j, step;
    char fcd_n = 'N';
    char fcd_t = 'T';
    int ione = 1, n1;
    int nrow, ncol, nmax;
    int IA=1 , JA=1;
    int *desca;


    nmax = lcr[iprobe].num_states;
    nrow = pmo.mxllda_lead[iprobe-1];
    ncol = pmo.mxlocc_lead[iprobe-1];
    desca = &pmo.desc_lead[(iprobe-1) * DLEN];


    n1 = nrow * ncol;

    /* allocate matrix and initialization  */
    alpha.r = 1.0;
    alpha.i = 0.0;
    beta.r = 0.0;
    beta.i = 0.0;
    my_malloc_init( tau,  nrow * ncol * 2, doublecomplex );
    my_malloc_init( taut,  nrow * ncol * 2, doublecomplex );
    my_malloc_init( t11,  nrow * ncol, doublecomplex );
    my_malloc_init( t12,  nrow * ncol, doublecomplex );
    my_malloc_init( ch0,  nrow * ncol, doublecomplex );
    my_malloc_init( ch1,  nrow * ncol, doublecomplex );

    my_malloc_init( ipiv, nrow + pmo.mblock, int );

    for (i = 0; i < nrow * ncol; i++)
    {
        ch0[i].r = Ha_eV * H00[i] - eneR * S00[i];
        ch0[i].i = -eneI * S00[i];
        ch1[i].r = Ha_eV * H01[i] - eneR * S01[i];
        ch1[i].i = -eneI * S01[i];
    }


    my_print_sum (ch0, n1, 1);
    my_print_sum (ch1, n1, 2);

    /* contruct the transfer matrix         */
    for (i = 0; i <  nrow * ncol; i++)
    {
        t12[i].r = -ch0[i].r;
        t12[i].i = -ch0[i].i;

    }

    /* t11 = (ene-ch0)^-1  */

    pmo_unitary_matrix(t11, &pmo.desc_lead[(iprobe-1) * DLEN]);

    PZGESV (&nmax, &nmax, t12, &IA, &JA, desca, ipiv, t11, &IA, &JA, desca, &info);



    if (info != 0)
    {
        printf ("Stransfer_p.c: error in PZGESV with INFO = %d in pe %d\n", info, pct.thispe);
        fflush (NULL);
        MPI_Finalize ();
        exit (0);
    }
    my_print_sum (t11, n1, 3);
    fflush(NULL);

    /* initialize intermediate t-matrices  */

    PZGEMM (&fcd_n, &fcd_t, &nmax, &nmax, &nmax, &alpha, t11, &IA, &JA, desca,
            ch1, &IA, &JA,  desca,  &beta, tau, &IA, &JA, desca);
    my_print_sum (tau, n1, 4);

    PZGEMM (&fcd_n, &fcd_n, &nmax, &nmax, &nmax, &alpha, t11, &IA, &JA, desca,
            ch1, &IA, &JA,  desca,  &beta, taut, &IA, &JA, desca);

    my_print_sum (taut, n1, 5);

    my_free(ch0);
    my_free(ch1);

    my_malloc_init( tsum, n1, doublecomplex );
    my_malloc_init( tsumt, n1, doublecomplex );
    my_malloc_init( s1, n1, doublecomplex );
    my_malloc_init( s2, n1, doublecomplex );

    ZCOPY (&n1, tau, &ione, tot, &ione);
    ZCOPY (&n1, taut, &ione, tott, &ione);
    ZCOPY (&n1, taut, &ione, tsum, &ione);
    ZCOPY (&n1, tau, &ione, tsumt, &ione);

    /*  iterative loop till convergence is achieved  */
    for (step = 0; step < MAX_STEP; step++)
    {

        PZGEMM (&fcd_n, &fcd_n, &nmax, &nmax, &nmax, &alpha, tau, &IA, &JA, desca,
                taut, &IA, &JA,  desca,  &beta, t11, &IA, &JA, desca);
        PZGEMM (&fcd_n, &fcd_n, &nmax, &nmax, &nmax, &alpha, taut, &IA, &JA, desca,
                tau, &IA, &JA,  desca,  &beta, t12, &IA, &JA, desca);

        my_print_sum (t11, n1, 6);
        my_print_sum (t12, n1, 7);

        pmo_unitary_matrix(s1, desca);
        pmo_unitary_matrix(s2, desca);

        for (i = 0; i < n1; i++)
        {
            s1[i].r += -t11[i].r - t12[i].r;
            s1[i].i += -t11[i].i - t12[i].i;
        }

        PZGESV (&nmax, &nmax, s1, &IA, &JA, desca, ipiv, s2, &IA, &JA, desca, &info);
        my_print_sum (s2, n1, 8);
        if (info != 0)
        {
            printf ("Stransfer.c: error in ZGESV with INFO = %d \n", info);
            fflush (NULL);
            MPI_Finalize ();
            exit (0);
        }

        PZGEMM (&fcd_n, &fcd_n, &nmax, &nmax, &nmax, &alpha, tau, &IA, &JA, desca,
                tau, &IA, &JA,  desca,  &beta, t11, &IA, &JA, desca);
        PZGEMM (&fcd_n, &fcd_n, &nmax, &nmax, &nmax, &alpha, taut, &IA, &JA, desca,
                taut, &IA, &JA,  desca,  &beta, t12, &IA, &JA, desca);
        PZGEMM (&fcd_n, &fcd_n, &nmax, &nmax, &nmax, &alpha, s2, &IA, &JA, desca,
                t11, &IA, &JA,  desca,  &beta, &tau[n1], &IA, &JA, desca);
        PZGEMM (&fcd_n, &fcd_n, &nmax, &nmax, &nmax, &alpha, s2, &IA, &JA, desca,
                t12, &IA, &JA,  desca,  &beta, &taut[n1], &IA, &JA, desca);


        my_print_sum (t11, n1, 9);
        my_print_sum (t12, n1, 10);
        my_print_sum (&tau[n1], n1, 11);
        my_print_sum (&taut[n1], n1, 12);

        PZGEMM (&fcd_n, &fcd_n, &nmax, &nmax, &nmax, &alpha, tsum, &IA, &JA, desca,
                &tau[n1], &IA, &JA,  desca,  &beta, t11, &IA, &JA, desca);
        PZGEMM (&fcd_n, &fcd_n, &nmax, &nmax, &nmax, &alpha, tsum, &IA, &JA, desca,
                &taut[n1], &IA, &JA,  desca,  &beta, s1, &IA, &JA, desca);

        my_print_sum (t11, n1, 301);
        my_print_sum (s1, n1, 302);

        ZCOPY (&n1, t11, &ione, s2, &ione);
        ZAXPY (&n1, &alpha, tot, &ione, s2, &ione);

        my_print_sum (s2, n1, 303);

        ZCOPY (&n1, s2, &ione, tot, &ione);
        ZCOPY (&n1, s1, &ione, tsum, &ione);

        PZGEMM (&fcd_n, &fcd_n, &nmax, &nmax, &nmax, &alpha, tsumt, &IA, &JA, desca,
                &taut[n1], &IA, &JA,  desca,  &beta, t11, &IA, &JA, desca);
        PZGEMM (&fcd_n, &fcd_n, &nmax, &nmax, &nmax, &alpha, tsumt, &IA, &JA, desca,
                &tau[n1], &IA, &JA,  desca,  &beta, s1, &IA, &JA, desca);

        my_print_sum (t11, n1, 304);
        my_print_sum (s1, n1, 305);

        ZCOPY (&n1, t11, &ione, s2, &ione);
        ZAXPY (&n1, &alpha, tott, &ione, s2, &ione);

        my_print_sum (s2, n1, 306);

        ZCOPY (&n1, s2, &ione, tott, &ione);
        ZCOPY (&n1, s1, &ione, tsumt, &ione);
        ZCOPY (&n1, &tau[n1], &ione, tau, &ione);
        ZCOPY (&n1, &taut[n1], &ione, taut, &ione);

        my_print_sum (tot, n1, 201);
        my_print_sum (tott, n1, 202);
        converge1 = 0.0;
        converge2 = 0.0;

        for (i = 0; i < n1; i++)
        {
            converge1 += sqrt (tau[ n1 + i].r * tau[n1 +i].r +
                    tau[ n1 + i].i * tau[n1 +i].i);
            converge2 += sqrt (taut[ n1 + i].r * taut[n1 +i].r +
                    taut[ n1 + i].i * taut[n1 +i].i);
        }

       comm_sums(&converge1, &ione, COMM_EN2);
       comm_sums(&converge2, &ione, COMM_EN2);

        if (converge1 < 1.0E-7 && converge2 < 1.0E-7)
            break;
    }

    if (converge1 > 1.0E-7 || converge2 > 1.0E-7)
    {
        printf ("bad t-matrix convergence\n");
        fflush (NULL);
        MPI_Finalize ();
        exit (0);
    }
    my_free(tau);
    my_free(taut);
    my_free(tsum);
    my_free(tsumt);
    my_free(t11);
    my_free(t12);
    my_free(s1);
    my_free(s2);
    my_free(ipiv);

}


static void my_print_sum (doublecomplex * matrix, int nmax, int index)
{
    int i;
    REAL a;
    a = 0.0;
    for (i = 0; i < nmax; i++)
        a += matrix[i].r + matrix[i].i;

    /*if(pct.thispe == 0) printf("\n sum of %f %d\n", a, index);  */
}
