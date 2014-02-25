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
#include <complex.h>

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"


void Sgreen_p (complex double * tot, complex double * tott, complex double *ch0, 
             complex double *ch1, complex double *g, int iprobe)
{

    complex double alpha, beta;
    int info;
    int *ipiv;
    int i;
    char fcd_n = 'N', fcd_t = 'T';
    int nrow, ncol, *desca, nmax, n1;
    int ione = 1;


    /* allocate matrix and initialization  */
    alpha = 1.0;
    beta = 1.0;

    nmax = lcr[iprobe].num_states;
    nrow = pmo.mxllda_lead[iprobe-1];
    ncol = pmo.mxlocc_lead[iprobe-1];
    n1 = nrow * ncol;
    desca = &pmo.desc_lead[ (iprobe-1) * DLEN];


    my_malloc_init( ipiv, nrow + pmo.mblock, int );


    PZGEMM (&fcd_n, &fcd_n, &nmax, &nmax, &nmax, &alpha, ch1, &ione, &ione, desca,
            tot, &ione, &ione, desca, &beta, ch0, &ione, &ione, desca);

    pmo_unitary_matrix(g, desca);


    PZGESV (&nmax, &nmax, ch0, &ione, &ione, desca, ipiv, g, &ione, &ione, desca, &info);
    if (info != 0)
    {
        printf ("Sgreen_p.c: error in ZGESV with INFO = %d \n", info);
        fflush (NULL);
        MPI_Finalize ();
        exit (0);
    }


    my_free( ipiv );

}
