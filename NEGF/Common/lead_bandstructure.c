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


char *get_num (char *str);
char *get_line (char *buf, FILE * fh);

void  pmo_unitary_matrix_double(double *, int *);

void PZHEGVX (int *, char *, char *, char *, int*, doublecomplex *, int*, int*,
        int *, doublecomplex *, int*, int*, int*, double *, double*, int*, int*,
        double *, int*, int*, double*,double*, doublecomplex *, int*,int*, int*,
        doublecomplex *, int*, double*, int*, int*, int*, 
        int*, int*, double*, int*);

void lead_bandstructure ()
{
    REAL de, emin, emax, E_image, kbt, current;


    int nmax, nC, nL, nR;
    int i, j, idx, idx1, idx2, idx3, E_POINTS;
    char jobz = 'N', uplo = 'U';
    FILE *file, *fhand;
    int *desca, mxllda, mxlocc, ndim;
    char tbuf[200], *tptr;

    doublecomplex *matH, *matS, *work, *z_vec;
    double *ener_band, *rwork, *eig_val;
    double kvec, coskvec, sinkvec, kmin, kmax, dk;
    double *temp, one = 1.0, zero = 0.0, *unitary_matrix;

    int itype = 1,ione = 1;
    char range = 'A';
    double VL = -100.0;
    double VU = 100.0;
    int IL = 1,izero = 0;
    int IU, nL1, nL2 ;       
    double tol = 1.0e-15, orfac = -1.0;
    doublecomplex *WORK, WORK_tmp ;
    int LWORK, LRWORK, *IWORK, IWORK_tmp, LIWORK, *IFAIL, *ICLUSTR;
    double *RWORK, RWORK_tmp, *GAP;
    int NB = pmo.mblock, NP0;

    int st1, st2, ik, kpoints, info;

    int ntot, iprobe;
    /* Open the input file for reading */

    read_cond_input (&emin, &emax, &E_POINTS, &E_image, &kbt, &kpoints);


    if (pct.thispe == 0)
    {
        printf ("\n band struture calculations for left lead\n\n ");
        printf ("lcr[1].num_states = %d \n", lcr[1].num_states);
        printf ("lcr[2].num_states = %d \n", lcr[2].num_states);
        printf (" kpoints= %d \n", kpoints);

    }


    desca = &pmo.desc_lead[0];
    nL = lcr[1].num_states;
    IU = nL;
    mxllda = pmo.mxllda_lead[0];
    mxlocc = pmo.mxlocc_lead[0];
    ndim = mxllda * mxlocc;


    ntot = 0;
    for (i = 0; i < ct.num_blocks; i++)
    {
        ntot += pmo.mxllda_cond[i] * pmo.mxlocc_cond[i];
    }
    for (i = 1; i < ct.num_blocks; i++)
    {
        ntot += pmo.mxllda_cond[i-1] * pmo.mxlocc_cond[i];
    }


    my_malloc_init( lcr[0].Htri, ntot, REAL );
    my_malloc_init( lcr[0].Stri, ntot, REAL );

    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx = pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];
        my_malloc_init( lcr[iprobe].H00, idx, REAL );
        my_malloc_init( lcr[iprobe].S00, idx, REAL );
        my_malloc_init( lcr[iprobe].H01, idx, REAL );
        my_malloc_init( lcr[iprobe].S01, idx, REAL );
    }

    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        i = cei.probe_in_block[iprobe - 1];
        idx = pmo.mxllda_cond[i] * pmo.mxlocc_lead[iprobe-1];
        my_malloc_init( lcr[iprobe].HCL, idx, REAL );
        my_malloc_init( lcr[iprobe].SCL, idx, REAL );
    }




    my_malloc_init( GAP, pmo.nrow * pmo.ncol, REAL );
    my_malloc_init( ICLUSTR, 2* pmo.nrow * pmo.ncol, int );
    my_malloc_init( IFAIL, nL, int );
    my_malloc_init( matH, mxllda * mxlocc, doublecomplex );
    my_malloc_init( matS, mxllda * mxlocc, doublecomplex );
    my_malloc_init( z_vec, mxllda * mxlocc, doublecomplex );
    my_malloc_init( unitary_matrix, mxllda * mxlocc, REAL );
    my_malloc_init( temp, mxllda * mxlocc, REAL );

    my_malloc_init( ener_band, kpoints * nL, REAL );

    my_malloc_init( eig_val, nL, REAL );

    pmo_unitary_matrix_double(unitary_matrix, desca);

    read_matrix_pp();

    kmin = 0.0;
    kmax = 4.0 * atan (1.0);
    dk = (kmax - kmin) / (kpoints - 1);

    NP0 = NUMROC( &nL, &NB, &izero, &izero, &pmo.nrow );
    LWORK = nL + max( NB * ( NP0 + 1 ), 3 );
    LRWORK = 9 * nL;
    LIWORK = 6 * nL;


    my_malloc_init( WORK, LWORK, doublecomplex );
    my_malloc_init( RWORK, LRWORK, double );
    my_malloc_init( IWORK, LIWORK, int );



    for (ik = pmo.myblacs; ik < kpoints; ik += pmo.npe_energy)
    {
        kvec = kmin + ik * dk;

        sinkvec = sin (kvec);
        coskvec = cos (kvec);

        /*  temp = trans (lcr[1].H01) */ 
        PDGEMM ("T", "N", &nL, &nL, &nL, &one, lcr[1].H01, &ione, &ione, desca,
                unitary_matrix, &ione, &ione, desca, &zero, temp, &ione, &ione, desca);

        for (st1 = 0; st1 < ndim; st1++)
        {

            matH[st1].r = lcr[1].H00[st1] + coskvec * (lcr[1].H01[st1] + temp[st1]);
            matH[st1].i = sinkvec * (lcr[1].H01[st1] - temp[st1]);

        }

        /*  temp = trans (lcr[1].S01) */ 
        PDGEMM ("T", "N", &nL, &nL, &nL, &one, lcr[1].S01, &ione, &ione, desca,
                unitary_matrix, &ione, &ione, desca, &zero, temp, &ione, &ione, desca);

        for (st1 = 0; st1 < ndim; st1++)
        {

            matS[st1].r = lcr[1].S00[st1] + coskvec * (lcr[1].S01[st1] + temp[st1]);
            matS[st1].i = sinkvec * (lcr[1].S01[st1] - temp[st1]);

        }


        PZHEGVX( &itype, &jobz, &range, &uplo, &nL, matH, &ione, &ione,
                desca, matS, &ione, &ione, desca, &VL, &VU, &IL, &IU,
                &tol, &nL1, &nL2, eig_val, &orfac,z_vec, &ione, &ione, desca,
                WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK,
                IFAIL, ICLUSTR, GAP, &info );

        if (info != 0)
        {

            printf ("\n info = %d %d %f \n", info, ik, kvec);
            exit (0);
        }

        for (st1 = 0; st1 < nL; st1++)
        {
            ener_band[ik * nL + st1] = 27.2116 * eig_val[st1];
        }


    }

    my_barrier ();
    idx = kpoints * nL;
    comm_sums (ener_band, &idx, COMM_EN1);

    my_barrier ();
    if (pct.thispe == 0)
    {
        file = fopen ("band.dat", "w");
        for (st1 = 0; st1 < nL; st1++)
        {

            for (ik = 0; ik < kpoints; ik++)
            {
                kvec = 1.0 * ik / (kpoints - 1);
                fprintf (file, " %15.8f,  %15.8f  \n", kvec, ener_band[ik * nL + st1]);
            }

            fprintf (file, "&\n");
        }
    }
    my_barrier ();

}
