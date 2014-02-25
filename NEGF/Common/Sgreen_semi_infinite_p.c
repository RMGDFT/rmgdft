/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
 *
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var_negf.h"
#include "LCR.h"
#include "pmo.h"


#define 	MAX_STEP 	40


void Sgreen_semi_infinite_p (complex double * green, complex double
        *ch00, complex double *ch01, int jprobe)
{

    double converge1, converge2, tem;
    complex double *chnn, *chtem;
    complex double alpha, beta, mone;
    int info;
    int *ipiv;
    int i, j, step;
    int ione = 1, n1;
    int maxrow, maxcol, *desca, nmax;
    char fcd1, fcd2;


    desca = &pmo.desc_lead[(jprobe -1) * DLEN];
    nmax = desca[2];

    fcd1 = 'N';
    fcd2 = 'T';


    maxrow = pmo.mxllda_lead[jprobe-1];
    maxcol = pmo.mxlocc_lead[jprobe-1];

    n1 = maxrow * maxcol;

    /* allocate matrix and initialization  */
    alpha = 1.0;
    beta = 0.0;
    mone = -1.0;

    my_malloc_init( chnn, n1, complex double );
    my_malloc_init( chtem, n1, complex double );
    my_malloc_init( ipiv, maxrow + pmo.mblock, int );


    /*  green = (e S00- H00)^-1  */

    ZCOPY (&n1, ch00, &ione, chnn, &ione);
    get_inverse_block_p (chnn, green, ipiv, desca);

    converge1 = 0.0;
    for (i = 0; i < n1; i++)
    {
        converge1 += cabs(green[i]) * cabs(green[i]);
    }

    comm_sums(&converge1, &ione, COMM_EN2);

    converge1 = sqrt (converge1);

    for (step = 0; step < MAX_STEP; step++)
    {

        /*  calculate chnn = ch00 - Hn+1, n * Gnn * Hn,n+1  */


        ZCOPY (&n1, ch00, &ione, chnn, &ione);
        PZGEMM (&fcd1, "N", &nmax, &nmax, &nmax, &alpha, ch01, &ione, &ione, desca,
                green, &ione, &ione, desca,  &beta, chtem, &ione, &ione, desca);
        PZGEMM ("N", &fcd2, &nmax, &nmax, &nmax, &mone, chtem, &ione, &ione, desca,
                ch01, &ione, &ione, desca, &alpha, chnn, &ione, &ione, desca);

        get_inverse_block_p (chnn, green, ipiv, desca);

        converge2 = 0.0;
        for (i = 0; i < n1; i++)
        {
            converge2 += cabs(green[i]) * cabs(green[i]) ;
            /* don't know what is the call to get norm */
        }
        comm_sums(&converge2, &ione, COMM_EN2);

        converge2 = sqrt (converge2);
        /* printf("\n  %d %f %f %16.8e converge \n", step, converge1, converge2, converge1-converge2); */

        tem = converge1 - converge2;
        tem = sqrt (tem * tem);

        if (tem < 1.0e-7)
            break;
        converge1 = converge2;
    }

    if (tem > 1.0e-7)
    {
        printf ("\n green not converge %f \n", tem);
        exit (0);
    }
    /*    printf("\n %d %f %f converge\n", step, eneR, eneI); */

    my_free( chnn );
    my_free( chtem );
    my_free( ipiv );
}
