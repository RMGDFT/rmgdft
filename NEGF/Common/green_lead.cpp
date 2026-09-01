#include "negf_prototypes.h"
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

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"
#include "GpuAlloc.h"
#include "transition.h"

#define 	MAX_STEP 	100


void green_lead (std::complex<double> *ch0_cpu, std::complex<double> *ch01_cpu, 
                std::complex<double> *ch10_cpu, std::complex<double> *green_cpu, int iprobe)
{

    double converge1, converge2;
    std::complex<double> *tau, *taut, *tsum, *tsumt, *t11, *t12, *s1;
    std::complex<double> one = 1.0, zero = 0.0, mone = -1.0;
    std::complex<double> *tot, *tott;
    std::complex<double> *temp_cpu,  *Imatrix_cpu;
    std::complex<double> *temp_gpu,  *Imatrix_gpu;
    std::complex<double> *temp_ptr,  *Imatrix_ptr;

    std::complex<double> *ch0_gpu, *ch01_gpu, *ch10_gpu, *green_gpu;
    std::complex<double> *ch0_ptr, *ch01_ptr, *ch10_ptr, *green_ptr;

    int step;
    int ione = 1, n1;
    int nrow, ncol, nmax;
    int *desca;

    nmax = lcr[iprobe].num_states;
    nrow = pmo.mxllda_lead[iprobe-1];
    ncol = pmo.mxlocc_lead[iprobe-1];
    desca = &pmo.desc_lead[(iprobe-1) * DLEN];


    n1 = nrow * ncol;

    /* allocate matrix and initialization  */
    size_t size =  n1 * sizeof(std::complex<double>);
    temp_cpu = (std::complex<double> *) RmgMallocHost( size * 11);
    Imatrix_cpu = (std::complex<double> *) RmgMallocHost( size );

    gpuMalloc((void **)&temp_gpu, size * 11);
    gpuMalloc((void **)&Imatrix_gpu, size );
    gpuMalloc((void **)&ch0_gpu, size );
    gpuMalloc((void **)&ch01_gpu, size );
    gpuMalloc((void **)&ch10_gpu, size );
    gpuMalloc((void **)&green_gpu, size );
    temp_ptr = MemoryPtrHostDevice(temp_cpu, temp_gpu);
    Imatrix_ptr = MemoryPtrHostDevice(Imatrix_cpu, Imatrix_gpu);
    ch0_ptr = MemoryPtrHostDevice(ch0_cpu, ch0_gpu);
    ch01_ptr = MemoryPtrHostDevice(ch01_cpu, ch01_gpu);
    ch10_ptr = MemoryPtrHostDevice(ch10_cpu, ch10_gpu);
    green_ptr = MemoryPtrHostDevice(green_cpu, green_gpu);

    tot  = &temp_ptr[0*n1];
    tott = &temp_ptr[1*n1];
    t11  = &temp_ptr[2*n1];
    t12  = &temp_ptr[3*n1];
    tsum = &temp_ptr[4*n1];
    tsumt= &temp_ptr[5*n1];
    s1   = &temp_ptr[6*n1];
    tau  = &temp_ptr[7*n1];
    taut = &temp_ptr[9*n1];
    

    pmo_unitary_matrix(Imatrix_cpu, desca);


    MemcpyHostDevice(size, Imatrix_cpu, Imatrix_gpu);
    MemcpyHostDevice(size, ch0_cpu, ch0_gpu);
    MemcpyHostDevice(size, ch01_cpu, ch01_gpu);
    MemcpyHostDevice(size, ch10_cpu, ch10_gpu);

    /* t11 = (ene-ch0)^-1  */

    zcopy_driver(n1, ch0_ptr, ione, t11, ione);
    matrix_inverse_driver(t11, desca);

    /* initialize intermediate t-matrices  */

    zgemm_driver ("N", "N", nmax, nmax, nmax, mone, t11, ione, ione, desca,
            ch10_ptr, ione, ione,  desca,  zero, tau, ione, ione, desca);

    zgemm_driver ("N", "N", nmax, nmax, nmax, mone, t11, ione, ione, desca,
            ch01_ptr, ione, ione,  desca,  zero, taut, ione, ione, desca);




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


        zcopy_driver (n1, Imatrix_ptr, ione, s1, ione);


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

        if(converge1 > 1.0E06 || converge2 > 1.0E06) 
            printf("\n WARNING Green function in green_lead.c diverging %d %e %e", step, converge1, converge2);
        if (converge1 < 1.0E-7 && converge2 < 1.0E-7)
            break;
    }

    if (converge1 > 1.0E-7 || converge2 > 1.0E-7)
    {
        rmg_printf ("bad t-matrix convergence\n");
        fflush (NULL);
        MPI_Finalize ();
        exit (0);
    }




    zgemm_driver ("N", "N", nmax, nmax, nmax, one, ch01_ptr, ione, ione, desca,
            tot, ione, ione, desca, one, ch0_ptr, ione, ione, desca);

    zcopy_driver(n1, ch0_ptr, ione, green_ptr, ione);
    matrix_inverse_driver(green_ptr, desca);

    MemcpyDeviceHost(size, green_gpu, green_cpu);
    RmgFreeHost(temp_cpu);
    RmgFreeHost(Imatrix_cpu);

    gpuFree(temp_gpu);
    gpuFree(Imatrix_gpu);
    gpuFree(ch0_gpu);
    gpuFree(ch01_gpu);
    gpuFree(ch10_gpu);
    gpuFree(green_gpu);
}


