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
#include <complex.h>

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"

#define 	MAX_STEP 	100


void *memory_ptr_host_device(void *ptr_host, void *ptr_device);
void matrix_inverse_driver(double *, int *);
void green_lead (complex double *ch0_host, complex double *ch01_host, 
                complex double *ch10_host, complex double *green_host, int iprobe)
{

    double converge1, converge2;
    complex double *tau, *taut, *tsum, *tsumt, *t11, *t12, *s1;
    complex double one = 1.0, zero = 0.0, mone = -1.0;
    complex double *green, *tot, *tott, *ch0, *ch01, *ch10; 
    complex double *temp_host, *Imatrix_host, *Imatrix;
    int info;
    int i, j, step;
    int ione = 1, n1;
    int nrow, ncol, nmax;
    int *desca;

    nmax = lcr[iprobe].num_states;
    nrow = pmo.mxllda_lead[iprobe-1];
    ncol = pmo.mxlocc_lead[iprobe-1];
    desca = &pmo.desc_lead[(iprobe-1) * DLEN];


    n1 = nrow * ncol;

    /* allocate matrix and initialization  */
    my_malloc_init( temp_host, n1 * 11, complex double );
    my_malloc_init( Imatrix_host, n1, complex double );

#if GPU_ENABLED
    if(pmo.ntot_low < 7 * n1)
    {
        printf("\n the central part is even smaller than lead, don't use gpu version \n");
        fflush(NULL);
        exit(0);
    }
#endif

    ch0 = memory_ptr_host_device(ch0_host, ct.gpu_Htri);
    ch10 = memory_ptr_host_device(ch10_host, &ct.gpu_Htri[n1]);
    ch01 = memory_ptr_host_device(ch01_host, ct.gpu_Hii);
    setvector_host_device (n1, sizeof(complex double), ch0_host, ione, ct.gpu_Htri, ione);
    setvector_host_device (n1, sizeof(complex double), ch10_host, ione, &ct.gpu_Htri[n1], ione);
    setvector_host_device (n1, sizeof(complex double), ch01_host, ione, ct.gpu_Hii, ione);

    tot  = memory_ptr_host_device(&temp_host[0*n1], &ct.gpu_Gtri[0*n1]);
    tott = memory_ptr_host_device(&temp_host[1*n1], &ct.gpu_Gtri[1*n1]);
    t11  = memory_ptr_host_device(&temp_host[2*n1], &ct.gpu_Gtri[2*n1]);
    t12  = memory_ptr_host_device(&temp_host[3*n1], &ct.gpu_Gtri[3*n1]);
    tsum = memory_ptr_host_device(&temp_host[4*n1], &ct.gpu_Gtri[4*n1]);
    tsumt= memory_ptr_host_device(&temp_host[5*n1], &ct.gpu_Gtri[5*n1]);
    s1   = memory_ptr_host_device(&temp_host[6*n1], &ct.gpu_Gtri[6*n1]);
    tau  = memory_ptr_host_device(&temp_host[7*n1], ct.gpu_Grow);
    taut = memory_ptr_host_device(&temp_host[9*n1], ct.gpu_Gcol);
    green= memory_ptr_host_device(green_host, ct.gpu_Gii);
    Imatrix = memory_ptr_host_device(Imatrix_host, ct.gpu_Imatrix);
    
    pmo_unitary_matrix(Imatrix_host, desca);
    setvector_host_device (n1, sizeof(complex double), Imatrix_host, ione, ct.gpu_Imatrix, ione);

    /* t11 = (ene-ch0)^-1  */

    zcopy_driver(n1, ch0, ione, t11, ione);
    matrix_inverse_driver(t11, desca);

    /* initialize intermediate t-matrices  */

    zgemm_driver ("N", "N", nmax, nmax, nmax, mone, t11, ione, ione, desca,
            ch10, ione, ione,  desca,  zero, tau, ione, ione, desca);

    zgemm_driver ("N", "N", nmax, nmax, nmax, mone, t11, ione, ione, desca,
            ch01, ione, ione,  desca,  zero, taut, ione, ione, desca);




    zcopy_driver (n1, tau, ione, tot, ione);
    zcopy_driver (n1, taut, ione, tott, ione);
    zcopy_driver (n1, taut, ione, tsum, ione);
    zcopy_driver (n1, tau, ione, tsumt, ione);

    /*  iterative loop till convergence is achieved  */
    for (step = 0; step < MAX_STEP; step++)
    {

        zgemm_driver ("N", "N", nmax, nmax, nmax, one, tau, ione, ione, desca,
                taut, ione, ione,  desca,  zero, t11, ione, ione, desca);
        zgemm_driver ("N", "N", nmax, nmax, nmax, one, taut, ione, ione, desca,
                tau, ione, ione,  desca,  zero, t12, ione, ione, desca);


        zcopy_driver (n1, Imatrix, ione, s1, ione);


        zaxpy_driver (n1, mone, t11, ione, s1, ione);
        zaxpy_driver (n1, mone, t12, ione, s1, ione);
        matrix_inverse_driver(s1, desca);

        zgemm_driver ("N", "N", nmax, nmax, nmax, one, tau, ione, ione, desca,
                tau, ione, ione,  desca,  zero, t11, ione, ione, desca);

        zgemm_driver ("N", "N", nmax, nmax, nmax, one, taut, ione, ione, desca,
                taut, ione, ione,  desca,  zero, t12, ione, ione, desca);

        zgemm_driver ("N", "N", nmax, nmax, nmax, one, s1, ione, ione, desca,
                t11, ione, ione,  desca,  zero, &tau[n1], ione, ione, desca);
        zgemm_driver ("N", "N", nmax, nmax, nmax, one, s1, ione, ione, desca,
                t12, ione, ione,  desca,  zero, &taut[n1], ione, ione, desca);



        zgemm_driver ("N", "N", nmax, nmax, nmax, one, tsum, ione, ione, desca,
                &tau[n1], ione, ione,  desca,  zero, t11, ione, ione, desca);
        zgemm_driver ("N", "N", nmax, nmax, nmax, one, tsum, ione, ione, desca,
                &taut[n1], ione, ione,  desca,  zero, s1, ione, ione, desca);

        zcopy_driver (n1, s1, ione, tsum, ione);


        zaxpy_driver (n1, one, t11, ione, tot, ione);


        zgemm_driver ("N", "N", nmax, nmax, nmax, one, tsumt, ione, ione, desca,
                &taut[n1], ione, ione,  desca,  zero, t11, ione, ione, desca);
        zgemm_driver ("N", "N", nmax, nmax, nmax, one, tsumt, ione, ione, desca,
                &tau[n1], ione, ione,  desca,  zero, s1, ione, ione, desca);
        zcopy_driver (n1, s1, ione, tsumt, ione);


        zaxpy_driver (n1, one, t11, ione, tott, ione);


        zcopy_driver (n1, &tau[n1], ione, tau, ione);
        zcopy_driver (n1, &taut[n1], ione, taut, ione);

        dzasum_driver(n1, &tau[n1], ione, &converge1);
        dzasum_driver(n1, &taut[n1], ione, &converge2);

        comm_sums(&converge1, &ione, COMM_EN2);
        comm_sums(&converge2, &ione, COMM_EN2);

        if(converge1 > 1.0E06 | converge2 > 1.0E06) 
            dprintf("\n WARNING Green function in green_lead.c diverging %e %e", converge1, converge2);
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




    zgemm_driver ("N", "N", nmax, nmax, nmax, one, ch01, ione, ione, desca,
            tot, ione, ione, desca, one, ch0, ione, ione, desca);

    zcopy_driver(n1, ch0, ione, green, ione);
    matrix_inverse_driver(green, desca);

    getvector_device_host (n1, sizeof(complex double), green, ione, green_host, ione);



    my_free(temp_host);
    my_free(Imatrix_host);

}


