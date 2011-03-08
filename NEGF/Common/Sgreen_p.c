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

#include "md.h"
#include "pmo.h"


void Sgreen_p (doublecomplex * tot, doublecomplex * tott, REAL * H00, REAL *H01,
             REAL * S00, REAL * S01, REAL eneR, REAL eneI, doublecomplex *g, 
             int iprobe)
{

    doublecomplex *ch0, *ch1 ;
    doublecomplex alpha, beta;
    int info;
    int *ipiv;
    int i;
    char fcd_n = 'N', fcd_t = 'T';
    int nrow, ncol, *desca, nmax, n1;
    int ione = 1;



    /* allocate matrix and initialization  */
    alpha.r = 1.0;
    alpha.i = 0.0;
    beta.r = 1.0;
    beta.i = 0.0;

    nmax = lcr[iprobe].num_states;
    nrow = pmo.mxllda_lead[iprobe-1];
    ncol = pmo.mxlocc_lead[iprobe-1];
    n1 = nrow * ncol;
    desca = &pmo.desc_lead[ (iprobe-1) * DLEN];


    my_malloc_init( ch0,  n1, doublecomplex);
    my_malloc_init( ch1,  n1, doublecomplex);

    my_malloc_init( ipiv, nrow + pmo.mblock, int );

    for (i = 0; i < n1; i++)
    {
        ch0[i].r = eneR * S00[i] - H00[i] * Ha_eV;
        ch0[i].i = eneI * S00[i];
        ch1[i].r = eneR * S01[i] - H01[i] * Ha_eV;
        ch1[i].i = eneI * S01[i];
    }


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


    my_free(ch0);
    my_free(ch1);
    my_free( ipiv );

}
