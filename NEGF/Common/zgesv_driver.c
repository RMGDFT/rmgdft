/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "make_conf.h"

#if GPU_ENABLED
//#include <magma.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <magma.h>
#endif

#include "my_scalapack.h"
#include "blas.h"
#include "RmgTimer.h"


// Ax = B

void zgesv_driver (complex double *A, int *desca,  complex double *B, int *descb )
{



    int d_ipiv, *ipiv, info, ione =1, izero = 0;
    int nprow, npcol, myrow, mycol;
    int ictxt = desca[1];
    int nn = desca[2], mb= desca[4], nb=desca[5];
    int nhrs = descb[3];
    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    void *RT5 = BeginRmgTimer("Ax=B");


    if(nprow*npcol <1) 
    {
        printf ("error in zgesv_driver nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }
#if GPU_ENABLED

    if(nprow*npcol != 1)
    {
        printf ("GPU ENALBED but nprow*npcol !=1  nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }
    
    d_ipiv = nn;


    ipiv = (int *) malloc(d_ipiv * sizeof(int));
    magma_zgesv_gpu (nn, nhrs, (magmaDoubleComplex *)A, nn, ipiv, (magmaDoubleComplex *)B, nn, &info);

    if (info != 0)
    {
        printf ("error in magma_zgesv with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }
    free(ipiv);
#else
    //  use scalapack if nprow * npcol > 1
    if(nprow*npcol > 1)  
    {

        d_ipiv = mb + numroc_( &nn, &mb, &myrow, &izero, &nprow);

        ipiv = (int *) malloc(d_ipiv * sizeof(int));
        pzgesv_ (&nn, &nhrs, (double *)A, &ione, &ione, desca, ipiv, (double *)B, &ione, &ione, descb, &info); 
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

    EndRmgTimer(RT5);

}


