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
#include "GpuAlloc.h"


#define 	MAX_STEP 	40

void Sgreen_semi_infinite_p (std::complex<double> *green, std::complex<double>
        *ch00, std::complex<double> *ch01, std::complex<double> *ch10, int jprobe)
{

    double converge1, converge2, tem;
    std::complex<double> *chtem;
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

    chtem = (std::complex<double> *)RmgMallocHost(n1 * sizeof(std::complex<double>));

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


    RmgFreeHost( chtem );
}
