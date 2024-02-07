#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if CUDA_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
    #include <cusolverDn.h>
    #if MAGMA_LIBS
        #include <magma.h>
    #endif
#endif


#include <complex>

#include "Scalapack.h"
#include "blas.h"
#include "RmgTimer.h"

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "Subdiag.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "Gpufuncs.h"
#include "blas.h"
#include "blacs.h"


#include "transition.h"
#include "RmgMatrix.h"



// Ax = B

void zgesv_driver (std::complex<double> *A, int *desca,  std::complex<double> *B, int *descb )
{



    int d_ipiv, *ipiv, info, ione =1, izero = 0;
    int nprow, npcol, myrow, mycol;
    int ictxt = desca[1];
    int nn = desca[2], mb= desca[4];
    int nhrs = descb[3];
    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    RmgTimer *RT5 = new RmgTimer("Ax=B");


    if(nprow*npcol <1) 
    {
        printf ("error in zgesv_driver nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }
#if CUDA_ENABLED || HIP_ENABLED

    if(nprow*npcol != 1)
    {
        printf ("GPU ENALBED but nprow*npcol !=1  nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }
    
    ZgetrftrsDriver(nn, nhrs, A, B);

#else
    //  use scalapack if nprow * npcol > 1
    if(nprow*npcol > 1)  
    {

        d_ipiv = mb + numroc_( &nn, &mb, &myrow, &izero, &nprow);

        ipiv = (int *) malloc(d_ipiv * sizeof(int));
        pzgesv (&nn, &nhrs, A, &ione, &ione, desca, ipiv, B, &ione, &ione, descb, &info); 
        if (info != 0)
        {
            printf ("error in pzgesv with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }

        free(ipiv);

    }
    else
    {

        d_ipiv = nn;
        ipiv = (int *) malloc(d_ipiv * sizeof(int));

        zgesv(&nn, &nhrs, (double *)A, &nn, ipiv, (double *)B, &nn, &info );
        if (info != 0)
        {
            printf ("error in zgesv with INFO = %d \n", info);
            fflush (NULL);
            exit (0);
        }

        free(ipiv);
    }
#endif

    delete RT5;

}

