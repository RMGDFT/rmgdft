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

void Sgreen_semi_infinite_p (std::complex<double> *green_cpu, std::complex<double>
        *ch00_cpu, std::complex<double> *ch01_cpu, std::complex<double> *ch10_cpu, int jprobe)
{

    double converge1, converge2, tem;
    std::complex<double> *chtem_cpu, *chtem_gpu, *chtem_ptr;
    std::complex<double> *ch00_gpu, *ch01_gpu, *ch10_gpu;
    std::complex<double> *ch00_ptr, *ch01_ptr, *ch10_ptr;
    std::complex<double> *green_gpu, *green_ptr;

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

    size_t size = n1 * sizeof(std::complex<double>);
    chtem_cpu = (std::complex<double> *)RmgMallocHost(size);

    gpuMalloc((void **)&chtem_gpu, size );
    gpuMalloc((void **)&ch00_gpu, size );
    gpuMalloc((void **)&ch01_gpu, size );
    gpuMalloc((void **)&ch10_gpu, size );
    gpuMalloc((void **)&green_gpu, size );
    chtem_ptr = MemoryPtrHostDevice(chtem_cpu, chtem_gpu);
    ch00_ptr = MemoryPtrHostDevice(ch00_cpu, ch00_gpu);
    ch01_ptr = MemoryPtrHostDevice(ch01_cpu, ch01_gpu);
    ch10_ptr = MemoryPtrHostDevice(ch10_cpu, ch10_gpu);
    green_ptr = MemoryPtrHostDevice(green_cpu, green_gpu);


    MemcpyHostDevice(size, ch00_cpu, ch00_gpu);
    MemcpyHostDevice(size, ch10_cpu, ch10_gpu);
    MemcpyHostDevice(size, ch01_cpu, ch01_gpu);



    /*  green = (e S00- H00)^-1  */

    zcopy_driver (n1, ch00_ptr, ione, green_ptr, ione);
    matrix_inverse_driver(green_ptr, desca);

    dzasum_driver(n1, green_ptr, ione, &converge1);

    comm_sums(&converge1, &ione, COMM_EN2);


    for (step = 0; step < MAX_STEP; step++)
    {

        /*  calculate chnn = ch00 - Hn+1, n * Gnn * Hn,n+1  */

        zgemm_driver ("N", "N", nmax, nmax, nmax, one, ch01_ptr, ione, ione, desca,
                green_ptr, ione, ione, desca,  zero, chtem_ptr, ione, ione, desca);
        zcopy_driver (n1, ch00_ptr, ione, green_ptr, ione);
        zgemm_driver ("N", "N", nmax, nmax, nmax, mone, chtem_ptr, ione, ione, desca,
                ch10_ptr, ione, ione, desca, one, green_ptr, ione, ione, desca);

        matrix_inverse_driver(green_ptr, desca);
        dzasum_driver(n1, green_ptr, ione, &converge2);

        comm_sums(&converge2, &ione, COMM_EN2);

        /* rmg_printf("\n  %d %f %f %16.8e converge \n", step, converge1, converge2, converge1-converge2); */

        tem = converge1 - converge2;
        tem = sqrt (tem * tem);

        if (tem < 1.0e-7)
            break;
        converge1 = converge2;
    }

    if (tem > 1.0e-7)
    {
        rmg_printf ("\n green not converge %f \n", tem);
        exit (0);
    }
    /*    rmg_printf("\n %d %f %f converge\n", step, eneR, eneI); */

    MemcpyDeviceHost(size, green_gpu, green_cpu);

    RmgFreeHost( chtem_cpu );
    gpuFree(chtem_gpu);
    gpuFree(ch00_gpu);
    gpuFree(ch10_gpu);
    gpuFree(ch01_gpu);
    gpuFree(green_gpu);

}
