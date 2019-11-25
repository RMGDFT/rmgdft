#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"

void *memory_ptr_host_device(void *ptr_host, void *ptr_device);
void matrix_inverse_driver(double *, int *);

// for GPU run, ct.gpu_Grow and ct.gpu_Gcol are already in device after Sgreen_c_noneq.c 
// Sgreen_c_noneq.c calls matrix_inverse_anyprobe.c with GPU 

void rho_munu (std::complex<double> * rho_mn_host, std::complex<double> * G_row_host, 
        std::complex<double> *G_col_host, std::complex<double> * gamma_host, int iprobe)
{

    std::complex<double> one = 1.0, zero=0.0, half=0.5;
    std::complex<double> *temp_host;
    std::complex<double> *Grow, *Gcol, *gamma, *temp, *rho_mn;
    int *desca, *descb, *descc, *descd, *desce, *descl;

    int i, n_green, n1, n2;

    int N, *ni, nL, N1;
    int maxrow, maxcol, ione =1;

    N = ct.num_blocks;
    ni = ct.block_dim;


    N1 = cei.probe_in_block[iprobe - 1];
    nL = ct.block_dim[N1];


    maxrow = 0;
    maxcol = 0;
    for (i = 0; i < N; i++)
    {
        maxrow = rmg_max(maxrow, pmo.mxllda_cond[i]);
        maxcol = rmg_max(maxcol, pmo.mxlocc_cond[i]);
    }

    n1 = maxrow * maxcol;
    my_malloc_init( temp_host, n1, std::complex<double> );

    gamma = memory_ptr_host_device(gamma_host, ct.gpu_Gii);
    temp = memory_ptr_host_device(temp_host, ct.gpu_temp);

    Grow = memory_ptr_host_device(G_row_host, ct.gpu_Grow);
    Gcol = memory_ptr_host_device(G_col_host, ct.gpu_Gcol);
    rho_mn = memory_ptr_host_device(rho_mn_host, ct.gpu_Gtri);


    n1 = nL * nL;
    setvector_host_device (n1, sizeof(std::complex<double>), gamma_host, ione, ct.gpu_Gii, ione);



    n_green = 0;
    descl = &pmo.desc_cond[ (N1  + N1 * ct.num_blocks ) * DLEN ];

    for (i = 0; i < N - 1; i++)
    {
        n1 = ni[i];
        n2 = ni[i + 1];
        desca = &pmo.desc_cond[ (i   +  N1   * ct.num_blocks ) * DLEN ];
        descb = &pmo.desc_cond[ (i   +  i    * ct.num_blocks ) * DLEN ];
        descc = &pmo.desc_cond[ (i+1 +  N1   * ct.num_blocks ) * DLEN ];
        descd = &pmo.desc_cond[ (i   + (i+1) * ct.num_blocks ) * DLEN ];

        /*  temp = G_i0 * Gamma  */
        zgemm_driver("N", "N", n1, nL, nL, one, &Grow[n_green], ione, ione, desca,
                gamma, ione, ione, descl, zero, temp, ione, ione, desca);

        /* rho_mn (i,i) = temp * G_i0^, the block (i,i) */

        zgemm_driver("N", "C", n1, n1, nL, one, temp, ione, ione, desca,
                &Grow[n_green], ione, ione, desca, zero, &rho_mn[pmo.diag_begin[i]], ione, ione, descb);

        /* rho_mn (i,i+1) = temp * G_i+10^, the block (i,i) */
        n_green += pmo.mxllda_cond[i] * maxcol;


        zgemm_driver("N", "C", n1, n2, nL, one, temp, ione, ione, desca,
                &Grow[n_green], ione, ione, descc, zero, &rho_mn[pmo.offdiag_begin[i]], ione, ione, descd);


    }

    /* calculate the last block  */
    n1 = ni[N - 1];

    desca = &pmo.desc_cond[ (N-1   +  N1   * ct.num_blocks ) * DLEN ];
    descb = &pmo.desc_cond[ (N-1   +  (N-1)   * ct.num_blocks ) * DLEN ];
    zgemm_driver("N", "N", n1, nL, nL, one, &Grow[n_green], ione, ione, desca,
            gamma, ione, ione, descl, zero, temp, ione, ione, desca);
    zgemm_driver("N", "C", n1, n1, nL, one, temp, ione, ione, desca,
            &Grow[n_green], ione, ione, desca, zero, &rho_mn[pmo.diag_begin[N-1]], ione, ione, descb);


    if(!ct.is_gamma)
    {

        n_green = 0;
        for (i = 0; i < N - 1; i++)
        {
            n1 = ni[i];
            n2 = ni[i + 1];
            desca = &pmo.desc_cond[ (N1 +  i   * ct.num_blocks ) * DLEN ];
            descb = &pmo.desc_cond[ (i  +  N1  * ct.num_blocks ) * DLEN ];
            descc = &pmo.desc_cond[ (i  +  i   * ct.num_blocks ) * DLEN ];
            descd = &pmo.desc_cond[ (N1 +  (i+1) * ct.num_blocks ) * DLEN ];
            desce = &pmo.desc_cond[ (i  +  (i+1) * ct.num_blocks ) * DLEN ];


            /*  temp = G_i0 * Gamma  */
            zgemm_driver("C", "N", n1, nL, nL, one, &Gcol[n_green], ione, ione, desca,
                    gamma, ione, ione, descl, zero, temp, ione, ione, descb);

            /* rho_mn (i,i) = temp * G_i0^, the block (i,i) */

            zgemm_driver("N", "N", n1, n1, nL, half, temp, ione, ione, descb,
                    &Gcol[n_green], ione, ione, desca, half, &rho_mn[pmo.diag_begin[i]], ione, ione, descc);

            /* rho_mn (i,i+1) = temp * G_i+10^, the block (i,i) */
            n_green += pmo.mxllda_cond[i] * maxcol;


            zgemm_driver("N", "N", n1, n2, nL, half, temp, ione,ione, descb,
                    &Gcol[n_green], ione, ione, descd, half, &rho_mn[pmo.offdiag_begin[i]], ione, ione, desce);


        }

        /* calculate the last block  */
        n1 = ni[N - 1];
        desca = &pmo.desc_cond[ (N1 +  (N-1) * ct.num_blocks ) * DLEN ];
        descb = &pmo.desc_cond[ (N-1 +  N1   * ct.num_blocks ) * DLEN ];
        descc = &pmo.desc_cond[ (N-1 + (N-1) * ct.num_blocks ) * DLEN ];

        zgemm_driver("C", "N", n1, nL, nL, one, &Gcol[n_green], ione, ione, desca,
                gamma, ione, ione, descl, zero, temp, ione, ione, descb);
        zgemm_driver("N", "N", n1, n1, nL, half, temp, ione, ione, descb,
                &Gcol[n_green], ione, ione, desca, half, &rho_mn[pmo.diag_begin[N-1]], ione, ione, descc);



    }

    getvector_device_host(pmo.ntot, sizeof( std::complex<double> ), ct.gpu_Gtri, ione, rho_mn_host, ione );


    int up_and_low = 0;
    green_kpoint_phase(rho_mn_host, ct.kp[pct.kstart].kpt[1], ct.kp[pct.kstart].kpt[2], up_and_low);

    my_free(temp_host);


}
