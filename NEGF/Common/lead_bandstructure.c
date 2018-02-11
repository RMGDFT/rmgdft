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
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"


char *get_num (char *str);
char *get_line (char *buf, FILE * fh);

void  pmo_unitary_matrix_double(double *, int *);

//void PZHEGVX (int *, char *, char *, char *, int*, complex double *, int*, int*,
//        int *, complex double *, int*, int*, int*, double *, double*, int*, int*,
//        double *, int*, int*, double*,double*, complex double *, int*,int*, int*,
//        complex double *, int*, double*, int*, int*, int*, 
//        int*, int*, double*, int*);

void lead_bandstructure ()
{
    double de, emin, emax, E_image, kbt, current;


    int nmax, nC, nL, nR;
    int i, j, idx, idx1, idx2, idx3, E_POINTS;
    char jobz = 'N', uplo = 'U';
    FILE *file, *fhand;
    int *desca, mxllda, mxlocc, ndim;
    char tbuf[200], *tptr;

    complex double *matH, *matS, *work, *z_vec;
    double *ener_band, *rwork, *eig_val;
    double kvec, coskvec, sinkvec, kmin, kmax, dk;
    double *temp,  *unitary_matrix;
    complex double one = 1.0, zero = 0.0;

    int itype = 1,ione = 1;
    char range = 'A';
    double VL = -100.0;
    double VU = 100.0;
    int IL = 1,izero = 0;
    int IU, nL1, nL2 ;       
    double tol = 1.0e-15, orfac = 1.0e-3;
    complex double *WORK, WORK_tmp ;
    int LWORK, LRWORK, *IWORK, IWORK_tmp, LIWORK, *IFAIL, *ICLUSTR;
    double *RWORK, RWORK_tmp, *GAP;
    int NB = pmo.mblock, NP0;

    int st1, st2, ik, kpoints[3], info;

    int ntot, iprobe;

    complex double *H10, *S10, *H01, *H00, *S01, *S00, *HCL, *SCL;

    complex double ctem1, ctem2;
    /* Open the input file for reading */

    read_cond_input (&emin, &emax, &E_POINTS, &E_image, &kbt, kpoints);


    if (pct.gridpe == 0)
    {
        printf ("\n band struture calculations for left lead\n\n ");
        printf ("lcr[1].num_states = %d \n", lcr[1].num_states);
        printf ("lcr[2].num_states = %d \n", lcr[2].num_states);
        printf (" kpoints= %d %d %d\n", kpoints[0], kpoints[1], kpoints[2]);

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

     int num_offdiag_yz = 9;

    my_malloc_init( lcr[0].Htri, ntot, double );
    my_malloc_init( lcr[0].Stri, ntot, double );

    my_malloc_init( lcr[0].Htri_yz, num_offdiag_yz *ntot, double );
    my_malloc_init( lcr[0].Stri_yz, num_offdiag_yz *ntot, double );



    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx = pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];
        my_malloc_init( lcr[iprobe].H00, idx, double );
        my_malloc_init( lcr[iprobe].S00, idx, double );
        my_malloc_init( lcr[iprobe].H01, idx, double );
        my_malloc_init( lcr[iprobe].S01, idx, double );


        my_malloc_init( lcr[iprobe].H00_yz, num_offdiag_yz * idx, double);
        my_malloc_init( lcr[iprobe].S00_yz, num_offdiag_yz * idx, double);
        my_malloc_init( lcr[iprobe].H01_yz, num_offdiag_yz * idx, double);
        my_malloc_init( lcr[iprobe].S01_yz, num_offdiag_yz * idx, double);


    }

    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        i = cei.probe_in_block[iprobe - 1];
        idx = pmo.mxllda_cond[i] * pmo.mxlocc_lead[iprobe-1];
        my_malloc_init( lcr[iprobe].HCL, idx, double );
        my_malloc_init( lcr[iprobe].SCL, idx, double );


        my_malloc_init( lcr[iprobe].HCL_yz, num_offdiag_yz *idx, double );
        my_malloc_init( lcr[iprobe].SCL_yz, num_offdiag_yz *idx, double);

    }




    my_malloc_init( GAP, pmo.nrow * pmo.ncol, double );
    my_malloc_init( ICLUSTR, 2* pmo.nrow * pmo.ncol, int );
    my_malloc_init( IFAIL, nL, int );
    my_malloc_init( matH, mxllda * mxlocc, complex double );
    my_malloc_init( matS, mxllda * mxlocc, complex double );
    my_malloc_init( z_vec, mxllda * mxlocc, complex double );
    my_malloc_init( unitary_matrix, mxllda * mxlocc, double );
    my_malloc_init( temp, mxllda * mxlocc, double );


    idx = mxllda * mxlocc;
    my_malloc_init( H00,  idx, complex double );
    my_malloc_init( H01,  idx, complex double );
    my_malloc_init( H10,  idx, complex double );
    my_malloc_init( S00,  idx, complex double );
    my_malloc_init( S01,  idx, complex double );
    my_malloc_init( S10,  idx, complex double );
    my_malloc_init( SCL,  idx, complex double );
    my_malloc_init( HCL,  idx, complex double );


    my_malloc_init( ener_band, kpoints[0] * nL, double );

    my_malloc_init( eig_val, nL, double );

    pmo_unitary_matrix_double(unitary_matrix, desca);

    read_matrix_pp();



    /* for center part, the orbital index is just same as input*/
    split_matrix_center ();

    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
        split_matrix_lead (iprobe);



    kmin = 0.0;
    kmax = 4.0 * atan (1.0);
    dk = (kmax - kmin) / (kpoints[0] - 1);

    NP0 = numroc( &nL, &NB, &izero, &izero, &pmo.nrow );
    LWORK = nL + rmg_max( NB * ( NP0 + 1 ), 3 );
    LRWORK = 9 * nL;
    LIWORK = 6 * nL;


    my_malloc_init( WORK, LWORK, complex double );
    my_malloc_init( RWORK, LRWORK, double );
    my_malloc_init( IWORK, LIWORK, int );

   double kvecy, kvecz;
    kvecy = 0.0;
    kvecz = 0.275 * 2.0 * 3.1415926;
   // kvecz = 0.0;

    iprobe = 1;
    idx = pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];

    matrix_kpoint_lead(S00, H00, S01, H01, SCL, HCL,  kvecy, kvecz, iprobe);

    //matrix_kpoint(idx, S01, lcr[iprobe].S01_yz, kvecy, kvecz);
    //matrix_kpoint(idx, H01, lcr[iprobe].H01_yz, kvecy, kvecz);
    //matrix_kpoint(idx, S00, lcr[iprobe].S00_yz, kvecy, kvecz);
    //matrix_kpoint(idx, H00, lcr[iprobe].H00_yz, kvecy, kvecz);

    desca = &pmo.desc_lead[ (iprobe-1) * DLEN];

    int numst = lcr[iprobe].num_states;

    pztranc(&numst, &numst, &one, S01, &ione, &ione, desca,
            &zero, S10, &ione, &ione, desca);
    pztranc(&numst, &numst, &one, H01, &ione, &ione, desca,
            &zero, H10, &ione, &ione, desca);




    for (ik = pmo.myblacs; ik < kpoints[0]; ik += pmo.npe_energy)
    {
        kvec = kmin + ik * dk;

        ctem1 = cexp(I*kvec);
        ctem2 = cexp(-I*kvec);

        /*  temp = trans (lcr[1].H01) */ 
 //       PDGEMM ("T", "N", &nL, &nL, &nL, &one, lcr[1].H01, &ione, &ione, desca,
 //               unitary_matrix, &ione, &ione, desca, &zero, temp, &ione, &ione, desca);

        for (st1 = 0; st1 < ndim; st1++)
        {

            matH[st1] = H00[st1] + ctem1 * H01[st1] + ctem2 *H10[st1];

        }

        /*  temp = trans (lcr[1].S01) */ 
//        PDGEMM ("T", "N", &nL, &nL, &nL, &one, lcr[1].S01, &ione, &ione, desca,
//                unitary_matrix, &ione, &ione, desca, &zero, temp, &ione, &ione, desca);

        for (st1 = 0; st1 < ndim; st1++)
        {

            matS[st1] = S00[st1] + ctem1 * S01[st1] + ctem2 * S10[st1];

        }


        pzhegvx( &itype, &jobz, &range, &uplo, &nL, (double *)matH, &ione, &ione,
                desca, (double *)matS, &ione, &ione, desca, &VL, &VU, &IL, &IU,
                &tol, &nL1, &nL2, eig_val, &orfac,(double *)z_vec, &ione, &ione, desca,
                (double *)WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK,
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
    idx = kpoints[0] * nL;
    comm_sums (ener_band, &idx, COMM_EN1);

    my_barrier ();
    if (pct.gridpe == 0)
    {
        file = fopen ("band.dat", "w");
        for (st1 = 0; st1 < nL; st1++)
        {

            for (ik = 0; ik < kpoints[0]; ik++)
            {
                kvec = 1.0 * ik / (kpoints[0] - 1);
                fprintf (file, " %15.8f,  %15.8f  \n", kvec, ener_band[ik * nL + st1]);
            }

            fprintf (file, "&\n");
        }
    }
    my_barrier ();

}
