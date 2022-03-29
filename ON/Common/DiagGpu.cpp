/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "init_var.h"
#include "RmgTimer.h"
#include "common_prototypes1.h"
#include "Scalapack.h"
#include "prototypes_on.h"

#include "blas.h"
#include "blas_driver.h"
#include "GpuAlloc.h"
#include "RmgGemm.h"
#include "RmgMatrix.h"



void DiagGpu(STATE *states, int numst, double *Hij_glob, double *Sij_glob, 
        double *rho_matrix_glob, double *theta_glob)
{

    RmgTimer  *RT0 = new RmgTimer("3-DiagGpu");

    int numst2 = numst * numst;
    int ione = 1;
    double *eigs= new double[numst];
    double *H_gpu, *S_gpu, *eigs_gpu, *eigvector_gpu, *work_gpu;
    size_t size = numst * numst * sizeof(double);
    gpuMalloc((void **)&H_gpu, size);
    gpuMalloc((void **)&S_gpu, size);
    gpuMalloc((void **)&eigvector_gpu, size);
    gpuMalloc((void **)&eigs_gpu, numst * sizeof(double) );


    MemcpyHostDevice(size, Hij_glob, H_gpu);
    MemcpyHostDevice(size, Sij_glob, S_gpu);

    RmgTimer *RT1b = new RmgTimer("3-DiagGpu: (S^-1)H");

    // H_gpu = S_gpu^(-1) H_gpu
    MemcpyHostDevice(size, Sij_glob, S_gpu);
    DgetrftrsDriver(numst, numst, S_gpu, H_gpu);
    double t1 = 2.0;
    dscal_driver(numst2, t1, H_gpu, ione);
    MemcpyDeviceHost(size, H_gpu, theta_glob);

    delete(RT1b);

    if(ct.num_gpu_devices == 1 )
    {
        int lwork = 3 * numst * numst + 8 * numst;
        lwork = std::max(lwork, 128000);
        gpuMalloc((void **)&work_gpu, lwork*sizeof(double));
        RmgTimer *RT1 = new RmgTimer("3-DiagGpu: Dsygvd");
        MemcpyHostDevice(size, Sij_glob, S_gpu);
        MemcpyHostDevice(size, Hij_glob, eigvector_gpu);
        DsygvdDriver((double *)eigvector_gpu, (double *)S_gpu, eigs_gpu, work_gpu, lwork, numst, numst);
        MemcpyDeviceHost(numst*sizeof(double), eigs_gpu, eigs);
        gpuFree(work_gpu);
        delete(RT1);
    }
    else if(ct.num_gpu_devices >1)
    {
        RmgTimer *RT1 = new RmgTimer("3-DiagGpu: DsygvdMg");
        if(pct.local_rank == 0 ) 
        {
            DsygvdMgDriver((double *)Hij_glob, (double *)Sij_glob, eigs, numst);
        }

        int factor = 2;
        if(ct.is_gamma) factor = 1;
        MPI_Bcast(Hij_glob, factor * numst*numst, MPI_DOUBLE, 0, pct.local_comm);
        MPI_Bcast(eigs, numst, MPI_DOUBLE, 0, pct.local_comm);

        MemcpyHostDevice(size, Hij_glob, eigvector_gpu);

        delete(RT1);
    }
    else
    {
    }


    for (int st1 = 0; st1 < numst; st1++)
    {
        states[st1].eig[0] = eigs[st1];
    }

    delete [] eigs;


    dcopy_driver(numst2, eigvector_gpu, ione, S_gpu, ione);
    for(int st1 = 0; st1 <  numst; st1++)
    {
        double alpha = states[st1].occupation[0];
        dscal_driver(numst, alpha, &S_gpu[st1 * numst], ione);
    }


    RmgTimer *RT3 = new RmgTimer("3-DiagGpu: gemm ");
    double one = 1.0, zero = 0.0;
    RmgGemm("N", "T", numst, numst, numst, one, eigvector_gpu, numst, S_gpu, numst, zero, H_gpu, numst);
    MemcpyDeviceHost(size, H_gpu, rho_matrix_glob);
    delete(RT3);


    gpuFree(H_gpu);
    gpuFree(S_gpu);
    gpuFree(eigvector_gpu);
    gpuFree(eigs_gpu);

    delete(RT0);
}

