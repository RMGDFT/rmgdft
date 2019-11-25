#include "negf_prototypes.h"
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


#define 	MAX_STEP 	40

void *memory_ptr_host_device(void *ptr_host, void *ptr_device);
void matrix_inverse_driver(double *, int *);


void Sgreen_semi_infinite_p (std::complex<double> * green_host, std::complex<double>
        *ch00_host, std::complex<double> *ch01_host, std::complex<double> *ch10_host, int jprobe)
{

    double converge1, converge2, tem;
    std::complex<double> *chtem, *chtem_host, *green, *ch00, *ch01, *ch10;
    std::complex<double> one=1.0, zero=0.0, mone=-1.0;
    int step;
    int ione = 1, n1;
    int maxrow, maxcol, *desca, nmax;


    desca = &pmo.desc_lead[(jprobe -1) * DLEN];
    nmax = desca[2];


    maxrow = pmo.mxllda_lead[jprobe-1];
    maxcol = pmo.mxlocc_lead[jprobe-1];

    n1 = maxrow * maxcol;

    /* allocate matrix and initialization  */

    my_malloc_init( chtem_host, n1, std::complex<double> );


    ch00 = memory_ptr_host_device(ch00_host, ct.gpu_Htri);
    ch10 = memory_ptr_host_device(ch10_host, &ct.gpu_Htri[n1]);
    ch01 = memory_ptr_host_device(ch01_host, ct.gpu_Hii);
    setvector_host_device (n1, sizeof(std::complex<double>), ch00_host, ione, ct.gpu_Htri, ione);
    setvector_host_device (n1, sizeof(std::complex<double>), ch10_host, ione, &ct.gpu_Htri[n1], ione);
    setvector_host_device (n1, sizeof(std::complex<double>), ch01_host, ione, ct.gpu_Hii, ione);

    chtem  = memory_ptr_host_device(chtem_host, ct.gpu_Gtri);
    green= memory_ptr_host_device(green_host, ct.gpu_Gii);



    /*  green = (e S00- H00)^-1  */

    zcopy_driver (n1, ch00, ione, green, ione);
    matrix_inverse_driver(green, desca);

    dzasum_driver(n1, green, ione, &converge1);

    comm_sums(&converge1, &ione, COMM_EN2);


    for (step = 0; step < MAX_STEP; step++)
    {

        /*  calculate chnn = ch00 - Hn+1, n * Gnn * Hn,n+1  */

        zgemm_driver ("N", "N", nmax, nmax, nmax, one, ch01, ione, ione, desca,
                green, ione, ione, desca,  zero, chtem, ione, ione, desca);
        zcopy_driver (n1, ch00, ione, green, ione);
        zgemm_driver ("N", "N", nmax, nmax, nmax, mone, chtem, ione, ione, desca,
                ch10, ione, ione, desca, one, green, ione, ione, desca);

        matrix_inverse_driver(green, desca);
        dzasum_driver(n1, green, ione, &converge2);

        comm_sums(&converge2, &ione, COMM_EN2);

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

    
    getvector_device_host (n1, sizeof(std::complex<double>), green, ione, green_host, ione);

    my_free( chtem_host );
}
