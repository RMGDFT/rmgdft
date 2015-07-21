/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"
#include "my_scalapack.h"

void *memory_ptr_host_device(void *ptr_host, void *ptr_device);
void matrix_inverse_driver(double *, int *);
void matrix_inverse_p (complex double * H_tri_host, complex double * G_tri_host)
{
/*  Calculate the inverse of a semi-tridiagonal complex matrix
 *
 *    H00   H01   0    0   ...		0
 *    H10   H11   H12  0  ...		0	
 *    0     H21   H22  H23 ....		0
 *    .     .     .    .		.	
 *    .     .     .    .		.
 *    .     .     .    .     Hn-1,n-1	Hn-1,n
 *    0     0     0          Hn,n-1     Hnn
 *
 *   H_tri: store the input matrix in the order of H00, H01, H11, H12, .... Hnn
 *   Hi,i+1  = transpose (Hi+1, i)
 *
 *   G_tri: output inverse matrix in the same order as H_tri
 *   N:   number of blocks
 *   ni:  dimension of each block
 *   for example Hii is a ni[i] x ni[i] matrix
 *               Hi,i+1 is a ni[i] x ni[i+1] matrix 
 */

    int  i, j, n1, n2, n3, n4, n5, n6, n7, n8;
    complex double *Gii, *G_tem, *G_col, *temp, *Hlower, *Hupper, *Gii0;
    complex double *Gii_host, *G_tem_host, *G_col_host, *temp_host;
    complex double *H_tri, *G_tri;
    complex double half, mone, one, zero;
    int ione = 1, *ntem_begin, *ncol_begin;
    int ntot_row, ntot_col;
    int maxrow, maxcol;
    int *desca, *descb, *descc, *descd;
    int *desce, *descf, *descg, *desch;
    int *desci, *descj, *desck, *descl;
    int *ni, N;

    complex double tttt[1];

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
        maxrow = max(maxrow, pmo.mxllda_cond[i]);
        maxcol = max(maxcol, pmo.mxlocc_cond[i]);

    }


    my_malloc_init( ntem_begin, ct.num_blocks, int);
    my_malloc_init( ncol_begin, ct.num_blocks, int);
    size_t n_alloc;
    n_alloc = maxrow * maxcol * sizeof(complex double);
    Gii_host = (complex double *) malloc(n_alloc);
    temp_host = (complex double *) malloc(n_alloc);

    n_alloc = ntot_row * maxcol * sizeof(complex double);
    G_tem_host = (complex double *) malloc(n_alloc);

    n_alloc = ntot_col * maxrow * sizeof(complex double);
    G_col_host = (complex double *) malloc(n_alloc);

    Gii = memory_ptr_host_device(Gii_host, ct.gpu_Gii);
    temp = memory_ptr_host_device(temp_host, ct.gpu_temp);
    G_tem = memory_ptr_host_device(G_tem_host, ct.gpu_Grow);
    G_col = memory_ptr_host_device(G_col_host, ct.gpu_Gcol);
    H_tri = memory_ptr_host_device(H_tri_host, ct.gpu_Htri);
    G_tri = memory_ptr_host_device(G_tri_host, ct.gpu_Gtri);

    /*
     *  ntem_begin: starting address of one column of G_tem for each block
     */

    ntem_begin[0] = 0;
    ncol_begin[0] = 0;
    for (i = 1; i < ct.num_blocks; i++)
    {
        ntem_begin[i] = ntem_begin[i - 1] + pmo.mxllda_cond[i - 1] * maxcol;
        ncol_begin[i] = ncol_begin[i - 1] + pmo.mxlocc_cond[i - 1] * maxrow;
    }
    

    /*  calculate the inverse of the first block  */

    setvector_host_device (pmo.ntot_low, sizeof(complex double), H_tri_host, ione, ct.gpu_Htri, ione);

    n1 = pmo.mxllda_cond[0] * pmo.mxlocc_cond[0];
    desca = &pmo.desc_cond[0];

    zcopy_driver (n1, H_tri, ione, Gii, ione);
    matrix_inverse_driver(Gii, desca);
    zcopy_driver (n1, Gii, ione, G_tri, ione);

    /*  iterate to get one more block  */

    for (i = 0; i < N - 1; i++)
    {
        /* get the interaction  Hi,i+1  from input H_tri 
         * Hupper is a pointer only
         */
        Hupper = &H_tri[pmo.offdiag_begin[i] ];
        Hlower = &H_tri[pmo.lowoffdiag_begin[i] ];
        Gii0 = &G_tri[pmo.diag_begin[i] ];

        /* Hii now has the matrix Hi+1,i+1  */

        n1 = pmo.mxllda_cond[i+1] * pmo.mxlocc_cond[i + 1]; 
        zcopy_driver (n1, &H_tri[pmo.diag_begin[i + 1]], ione, Gii, ione);


        /* calculate Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1  */

        n1 = ni[i + 1];
        n2 = ni[i];

        desca = &pmo.desc_cond[ (i   + (i+1) * ct.num_blocks) * DLEN];
        descb = &pmo.desc_cond[ (i   +     i * ct.num_blocks) * DLEN];
        descc = &pmo.desc_cond[ (i+1 +     i * ct.num_blocks) * DLEN];
        descd = &pmo.desc_cond[ (i+1 + (i+1) * ct.num_blocks) * DLEN];

        zgemm_driver ("N", "N", n1, n2, n2, one, Hlower, ione, ione, descc,
                Gii0, ione, ione, descb, zero, temp, ione, ione, descc);
        zgemm_driver ("N", "N", n1, n1, n2, mone, temp, ione, ione, descc,
                Hupper, ione, ione, desca, one, Gii, ione, ione, descd);


        /* now Hii store the matrix Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1
         * Gi+1,i+1, stored in Gii, = Hii^(-1)
         */

        matrix_inverse_driver(Gii, descd);

        n1 = pmo.mxllda_cond[i + 1] * pmo.mxlocc_cond[i + 1];
        zcopy_driver (n1, Gii, ione, &G_tri[pmo.diag_begin[i + 1]], ione);

        /*  Gi,i+1 =  Gii^0 * Hi,i+1 * Gi+1,i+1  */

        n1 = ni[i];
        n2 = ni[i + 1];
        n3 = pmo.offdiag_begin[i];

        /* temp = Hi+1,i * Gii^0  eq. 11 for j = I term*/
        zgemm_driver ("N", "N", n2, n1, n1, mone, Hlower, ione, ione, descc,
                Gii0, ione, ione, descb, zero, temp, ione, ione, descc);

        zgemm_driver ("N", "N", n2, n1, n2, one, &G_tri[pmo.diag_begin[i+1]], ione, ione, descd,
                temp, ione, ione, descc, zero, &G_tri[pmo.lowoffdiag_begin[i]], ione, ione, descc);

        n4 = pmo.mxllda_cond[i+1] * pmo.mxlocc_cond[i];
        zcopy_driver(n4, &G_tri[pmo.lowoffdiag_begin[i]], ione, &G_col[ncol_begin[i]], ione); 

        /* temp = Gii^0 * Hi,i+1 */
        zgemm_driver ("N", "N", n1, n2, n1, mone, Gii0, ione, ione, descb,
                Hupper, ione, ione, desca, zero, temp, ione, ione, desca);

        /* G(i,i+1) = temp * G(i+1,i+1)  also == G_tem(i,i+1)  */
        zgemm_driver ("N", "N", n1, n2, n2, one, temp, ione, ione, desca,
                Gii, ione, ione, descd, zero, &G_tri[n3], ione, ione, desca);

        n4 = pmo.mxllda_cond[i] * pmo.mxlocc_cond[i + 1];
        zcopy_driver (n4, &G_tri[n3], ione, &G_tem[ntem_begin[i]], ione);

        /* update Gii  */

        zgemm_driver ("N", "N", n1, n1, n2, one, temp, ione, ione, desca,
                &G_tri[pmo.lowoffdiag_begin[i]], ione, ione, descc, one, &G_tri[pmo.diag_begin[i]], ione, ione, descb);

        for (j = i - 1; j >= 0; j--)
        {

            n1 = ni[j];         /* dimension of block j */
            n2 = ni[i];         /* dimension of block i */
            n3 = ni[i + 1];     /* dimension of block i+1 */
            n4 = ntem_begin[j]; /* starting address of Gtem(j,*) in G_tem */
            n5 = pmo.diag_begin[j];    /* starting address of G(j,j) in G_tri */
            n6 = pmo.offdiag_begin[j];    /* starting address of G(j,j+1) in G_tri */
            n7 = ntem_begin[j + 1];     /* starting address of Gtem(j+1,*) in G_tem */
            n8 = ni[j + 1];     /* dimension of block j+1 */

            desca = &pmo.desc_cond[ ( j   +  i    * ct.num_blocks) * DLEN];
            descb = &pmo.desc_cond[ ( i   + (i+1) * ct.num_blocks) * DLEN];
            descc = &pmo.desc_cond[ ( j   + (i+1) * ct.num_blocks) * DLEN];
            descd = &pmo.desc_cond[ ( i+1 +  i    * ct.num_blocks) * DLEN];
            desce = &pmo.desc_cond[ ( i   +  j    * ct.num_blocks) * DLEN];
            descf = &pmo.desc_cond[ ( i+1 +  j    * ct.num_blocks) * DLEN];
            descg = &pmo.desc_cond[ ( i+1 + (i+1) * ct.num_blocks) * DLEN];
            desch = &pmo.desc_cond[ ( j   +  j    * ct.num_blocks) * DLEN];
            desci = &pmo.desc_cond[ ( i+1 + (j+1) * ct.num_blocks) * DLEN];
            descj = &pmo.desc_cond[ ( j   + (j+1) * ct.num_blocks) * DLEN];
            desck = &pmo.desc_cond[ ( j+1 + (i+1) * ct.num_blocks) * DLEN];
            descl = &pmo.desc_cond[ ( j+1 +  j    * ct.num_blocks) * DLEN];


            /* gpu_Gii = -Hi+1,i * G0(i,j)   */

            zgemm_driver ("N", "N", n3, n1, n2, mone, Hlower, ione, ione, descd,
                    &G_col[ncol_begin[j]], ione, ione, desce, zero, Gii, ione, ione, descf);

            /*  G(I+1, j) = G(I+1, I+1) * gpu_Gii  eq. 11 */
            zgemm_driver ("N", "N", n3, n1, n3, one, &G_tri[pmo.diag_begin[i+1]], ione, ione, descg,
                    Gii, ione, ione, descf, zero, &G_col[ncol_begin[j]], ione, ione, descf);


            /* temp = -G0(j,i) * Hi,i+1  */
            zgemm_driver ("N", "N", n1, n3, n2, mone, &G_tem[n4], ione, ione, desca,
                    Hupper, ione, ione, descb, zero, temp, ione, ione, descc);

            /* G0(j, i+1) = temp * G(i+1,i+1) eq 8 */

            zgemm_driver ("N", "N", n1, n3, n3, one, temp, ione, ione, descc,
                    &G_tri[pmo.diag_begin[i+1]], ione, ione, descg, zero, &G_tem[n4], ione, ione, descc);

            /* G(j,j) = G0(j,j) + temp * G(i+1,j) eq.9  */

            zgemm_driver ("N", "N", n1, n1, n3, one, temp, ione, ione, descc, 
                    &G_col[ncol_begin[j]], ione, ione, descf, one, &G_tri[n5], ione, ione, desch);

            /* G(j,j+1) = G0(j,j+1) + temp * G(i+1,j+1) eq. 10 */
            zgemm_driver ("N", "N", n1, n8, n3, one, temp, ione, ione, descc, 
                    &G_col[ncol_begin[j+1]], ione, ione, desci, one, &G_tri[n6],
                    ione, ione, descj);

            /* G(j,j+1) = G0(j,j+1) + temp * G(i+1,j+1)  eq. 12 */
            zgemm_driver ("N", "N", n8, n1, n3, one, &G_tem[n7], ione, ione, desck, 
                    Gii, ione, ione, descf, one, &G_tri[pmo.lowoffdiag_begin[j]],
                    ione, ione, descl);

        }                       /* end for (j--) */

    }                           /* end  for(i = 0; i < N-1; i++) */

    getvector_device_host (pmo.ntot_low, sizeof(complex double),ct.gpu_Gtri,ione, G_tri_host, ione);

    int up_and_low = 1;
    green_kpoint_phase(G_tri_host, ct.kp[pct.kstart].kpt[1], ct.kp[pct.kstart].kpt[2], up_and_low);

    for(i = 0; i < ct.num_blocks - 1; i++)
    {
        n1 = ct.block_dim[i];
        n2 = ct.block_dim[i+1];
        n3 = pmo.offdiag_begin[i];
        n4 = pmo.lowoffdiag_begin[i];

        desca = &pmo.desc_cond[ ( i +  (i+1)  * ct.num_blocks) * DLEN];
        descc = &pmo.desc_cond[ ( i+1 +  i    * ct.num_blocks) * DLEN];

        pztranu_(&n1, &n2, &half, &G_tri_host[n4], &ione, &ione, descc, 
                &half, &G_tri_host[n3], &ione, &ione, desca);

    }



    my_free( ntem_begin );
    my_free( ncol_begin );
    free( Gii_host );
    free( temp_host );
    free( G_tem_host);
    free( G_col_host );
}
