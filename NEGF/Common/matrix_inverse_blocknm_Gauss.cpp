#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"
#include "GpuAlloc.h"



void matrix_inverse_blocknm_Gauss (std::complex<double> * H_tri_cpu, std::complex<double> *G_tri_cpu, 
        int m, int n, std::complex<double> * Green_C_cpu)
{
/*  Calculate the inverse of a semi-tridiagonal complex matrix
 *
 *    H00   H01   0    0   ...      0
 *    H10   H11   H12  0  ...       0   
 *    0     H21   H22  H23 ....     0
 *    .     .     .    .        .   
 *    .     .     .    .        .
 *    .     .     .    .     Hn-1,n-1   Hn-1,n
 *    0     0     0          Hn,n-1     Hnn
 *
 *   H_tri: store the input matrix in the order of H00, H01, H11, H12, .... Hnn
 *   Hi,i+1  = transpose (Hi+1, i)
 *
 *   Green_C: block (m2,m1) will be returned which will be used to
 *   calculate the transmission.
 *   N:   number of blocks
 *   ni:  dimension of each block
 *   for example Hii is a ni[i] x ni[i] matrix
 *               Hi,i+1 is a ni[i] x ni[i+1] matrix 
 */

    int  i, n1, n2;
    std::complex<double> *Hlower, *Hupper;
    std::complex<double> *Gdiag_cpu, *Gdiag_gpu, *Gdiag_ptr;
    std::complex<double> *Gii_cpu, *Gii_gpu, *Gii_ptr;
    std::complex<double> *H_tri_gpu, *H_tri_ptr;
    std::complex<double> *G_tri_gpu, *G_tri_ptr;
    std::complex<double> *Green_C_gpu, *Green_C_ptr;

    std::complex<double> half, mone, one, zero;
    int ione = 1;
    int *ndiag_begin;
    int ntot_row, ntot_col;
    int maxrow, maxcol;
    int *desca, *descb, *descc, *descd;
    int ncopy, N, *ni;

    ni = ct.block_dim;
    N = ct.num_blocks;
    mone = -1.0;
    one = 1.0;
    half = 0.5;
    zero = 0.0;

    /*  find the maximum dimension of the blocks  */


    ntot_row = 0;
    ntot_col = 0;
    maxrow=0;
    maxcol=0;
    for (i = 0; i < ct.num_blocks; i++)
    {
        ntot_row += pmo.mxllda_cond[i];
        ntot_col += pmo.mxlocc_cond[i];
        maxrow = rmg_max(maxrow, pmo.mxllda_cond[i]);
        maxcol = rmg_max(maxcol, pmo.mxlocc_cond[i]);

    }


    my_malloc_init( ndiag_begin, ct.num_blocks, int);
    size_t n_alloc;
    n_alloc = maxrow * maxcol * sizeof(std::complex<double>);
    Gii_cpu = Green_C_cpu;
    gpuMalloc((void **)&Gii_gpu, n_alloc );
    Gii_ptr = MemoryPtrHostDevice(Gii_cpu, Gii_gpu);

    gpuMalloc((void **)&H_tri_gpu, pmo.ntot_low * sizeof(std::complex<double>) );
    gpuMalloc((void **)&G_tri_gpu, pmo.ntot_low * sizeof(std::complex<double>) );
    H_tri_ptr = MemoryPtrHostDevice(H_tri_cpu, H_tri_gpu);
    G_tri_ptr = MemoryPtrHostDevice(G_tri_cpu, G_tri_gpu);

    MemcpyHostDevice(pmo.ntot_low * sizeof(std::complex<double>), H_tri_cpu, H_tri_gpu);


    n_alloc = 0;
    for (i = 0; i < ct.num_blocks; i++)
    {
        n_alloc += pmo.mxllda_cond[i] * pmo.mxlocc_cond[i];
    }

    Gdiag_cpu = (std::complex<double> *) RmgMallocHost(n_alloc * sizeof(std::complex<double>));
    gpuMalloc((void **)&Gdiag_gpu, n_alloc * sizeof(std::complex<double>) );
    Gdiag_ptr = MemoryPtrHostDevice(Gdiag_cpu, Gdiag_gpu);


    /*
     *  ndiag_begin[i]:  pointer address for i-th diagonal block in Gdiag
     */

    ndiag_begin[0] = 0;
    for (i = 1; i < ct.num_blocks; i++)
    {
        ndiag_begin[i] = ndiag_begin[i - 1] + pmo.mxlocc_cond[i - 1] * pmo.mxllda_cond[i - 1];
    }
    

    /*  calculate the inverse of the first block  */

    //  right side Gauss elimination  

    ncopy = pmo.mxllda_cond[N-1] * pmo.mxlocc_cond[N-1];

    zcopy_driver (ncopy, &H_tri_ptr[pmo.diag_begin[N-1]], ione, &Gdiag_ptr[ndiag_begin[N-1]], ione);


    for (i = N-1; i > m; i--)
    {
        /* get the interaction  Hi,i+1  from input H_tri 
         * Hupper is a pointer only  Hi-1, i
         * Hlower is a pointer only  Hi, i-1
         */
        Hupper = &H_tri_ptr[pmo.offdiag_begin[i-1] ];
        Hlower = &H_tri_ptr[pmo.lowoffdiag_begin[i-1] ];


        desca = &pmo.desc_cond[ (i   +     i * ct.num_blocks) * DLEN];
        descb = &pmo.desc_cond[ (i-1 +     i * ct.num_blocks) * DLEN];
        descc = &pmo.desc_cond[ (i   + (i-1) * ct.num_blocks) * DLEN];
        descd = &pmo.desc_cond[ (i-1 + (i-1) * ct.num_blocks) * DLEN];

        n1 = ni[i-1];
        n2 = ni[i];

        //Ci = (Dii)^-1 * Hi,i-1
        ncopy = pmo.mxllda_cond[i] * pmo.mxlocc_cond[i]; 
        zcopy_driver (ncopy, &Gdiag_ptr[ndiag_begin[i]], ione, Gii_ptr, ione);


        ncopy = pmo.mxllda_cond[i] * pmo.mxlocc_cond[i-1]; 
        zcopy_driver (ncopy, Hlower, ione, &G_tri_ptr[pmo.lowoffdiag_begin[i-1]], ione);
        zgesv_driver(Gii_ptr, desca, &G_tri_ptr[pmo.lowoffdiag_begin[i-1]], descc);


        //  Di+1, i+1 = Hi+1,i+1 +Ci * Hi,i+1

        ncopy = pmo.mxllda_cond[i-1] * pmo.mxlocc_cond[i - 1]; 
        zcopy_driver (ncopy, &H_tri_ptr[pmo.diag_begin[i - 1]], ione, &Gdiag_ptr[ndiag_begin[i-1]], ione);

        zgemm_driver ("N", "N", n1, n1, n2, mone, Hupper, ione, ione, descb, 
                &G_tri_ptr[pmo.lowoffdiag_begin[i-1]], ione, ione, descc,
                one, &Gdiag_ptr[ndiag_begin[i-1]], ione, ione, descd);
    }

     //  left side Gauss elimination  

    ncopy = pmo.mxllda_cond[0] * pmo.mxlocc_cond[0];

    zcopy_driver (ncopy, H_tri_ptr, ione, G_tri_ptr, ione);


    for (i = 0; i < m; i++)
    {
        /* get the interaction  Hi,i+1  from input H_tri 
         * Hupper is a pointer only  Hi, i+1
         * Hlower is a pointer only  Hi+1, i
         */
        Hupper = &H_tri_ptr[pmo.offdiag_begin[i] ];
        Hlower = &H_tri_ptr[pmo.lowoffdiag_begin[i] ];


        desca = &pmo.desc_cond[ (i   +     i * ct.num_blocks) * DLEN];
        descb = &pmo.desc_cond[ (i+1 +     i * ct.num_blocks) * DLEN];
        descc = &pmo.desc_cond[ (i   + (i+1) * ct.num_blocks) * DLEN];
        descd = &pmo.desc_cond[ (i+1 + (i+1) * ct.num_blocks) * DLEN];

        n1 = ni[i+1];
        n2 = ni[i];

        ncopy = pmo.mxllda_cond[i] * pmo.mxlocc_cond[i]; 
        zcopy_driver (ncopy, &G_tri_ptr[pmo.diag_begin[i]], ione, Gii_ptr, ione);

        ncopy = pmo.mxllda_cond[i] * pmo.mxlocc_cond[i+1]; 
        zcopy_driver (ncopy, Hupper, ione, &G_tri_ptr[pmo.offdiag_begin[i]], ione);
        //  Ci = -(Di,i)^-1 * Hi,i+1
        zgesv_driver (Gii_ptr, desca, &G_tri_ptr[pmo.offdiag_begin[i]], descc);

       //  Di+1, i+1 = Hi+1,i+1 +Ci * Hi,i+1

        ncopy = pmo.mxllda_cond[i+1] * pmo.mxlocc_cond[i + 1]; 
        zcopy_driver (ncopy, &H_tri_ptr[pmo.diag_begin[i + 1]], ione, &G_tri_ptr[pmo.diag_begin[i+1]], ione);

        zgemm_driver ("N", "N", n1, n1, n2, mone, Hlower, ione, ione, descb,
                &G_tri_ptr[pmo.offdiag_begin[i]], ione, ione, descc,
                one, &G_tri_ptr[pmo.diag_begin[i+1]], ione, ione, descd);
    }





    //  now G_tri  diag part store Dii_L from left
    //      Gdiag  store Dii_R from right
    //   G_tri  upper off diag store Ci_L from left
    //   G_tri  low off diag store Ci_R from right 

    //calculating  diag blocks of G_tri = (-Hii + Dii_L + Dii_R)^-1


    i = m;
    {
        n1 = ni[i];
        desca = &pmo.desc_cond[ (i   +     i * ct.num_blocks) * DLEN];

        ncopy = pmo.mxllda_cond[i] * pmo.mxlocc_cond[i]; 
        zaxpy_driver (ncopy, one, &Gdiag_ptr[ndiag_begin[i]], ione, &G_tri_ptr[pmo.diag_begin[i]], ione);
        zaxpy_driver (ncopy, mone, &H_tri_ptr[pmo.diag_begin[i]], ione, &G_tri_ptr[pmo.diag_begin[i]], ione);
        matrix_inverse_driver(&G_tri_ptr[pmo.diag_begin[i]], desca);

    }


    ncopy = pmo.mxllda_cond[m] * pmo.mxlocc_cond[m]; 
    zcopy_driver (ncopy, &G_tri_ptr[pmo.diag_begin[m]], ione, Gii_ptr, ione);

    int n0 = ni[m];


    // for case m M b
    for(i = m; i < n; i++)
    {
        n1 = ni[i];
        n2 = ni[i+1];
        desca = &pmo.desc_cond[ ((i+1) + i * ct.num_blocks) * DLEN];
        descb = &pmo.desc_cond[ ((i+1) + m * ct.num_blocks) * DLEN];
        descc = &pmo.desc_cond[ ( i    + m * ct.num_blocks) * DLEN];

        zgemm_driver ("N", "N", n2, n0, n1, mone, &G_tri_ptr[pmo.lowoffdiag_begin[i]], ione, ione, desca,
                Gii_ptr, ione, ione, descc, zero, Gdiag_ptr, ione, ione, descb);


        ncopy = pmo.mxllda_cond[i+1] * pmo.mxlocc_cond[m]; 
        zcopy_driver (ncopy, Gdiag_ptr, ione, Gii_ptr, ione);
    }

    for(i = m; i > n; i--)
    {
        n1 = ni[i];
        n2 = ni[i-1];
        desca = &pmo.desc_cond[ ((i-1) + i * ct.num_blocks) * DLEN];
        descb = &pmo.desc_cond[ ((i-1) + m * ct.num_blocks) * DLEN];
        descc = &pmo.desc_cond[ ( i    + m * ct.num_blocks) * DLEN];

        zgemm_driver ("N", "N", n2, n0, n1, mone, &G_tri_ptr[pmo.offdiag_begin[i-1]], ione, ione, desca,
                Gii_ptr, ione, ione, descc, zero, Gdiag_ptr, ione, ione, descb);


        ncopy = pmo.mxllda_cond[i-1] * pmo.mxlocc_cond[m]; 
        zcopy_driver (ncopy, Gdiag_ptr, ione, Gii_ptr, ione);
    }


    ncopy = pmo.mxllda_cond[n] * pmo.mxlocc_cond[m]; 
    MemcpyDeviceHost(ncopy * sizeof(std::complex<double>), Gii_gpu, Gii_cpu);


    my_free( ndiag_begin );
    RmgFreeHost( Gdiag_cpu );
    gpuFree(Gii_gpu);
    gpuFree(Gdiag_gpu);
    gpuFree(G_tri_gpu);
    gpuFree(H_tri_gpu);

}

