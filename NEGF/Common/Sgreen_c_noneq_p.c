/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
 *  Calculate the right-up block for transmission calculations.
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


void Sgreen_c_noneq_p (double *Htri, double *Stri, complex double * sigma,
                     int *sigma_idx, complex double ene, complex double *Green_C_row, 
                     complex double *Green_C_col, int nC, int iprobe)
{
/*   H00, S00: nC * nC real matrix
 *   sigma:  nC * nC complex matrix
 *   Green_C  nC * nC complex matrix , output green function
 *   nC: number of states in conductor region
 */

    complex double *H_tri, *H_whole, *H_inv;

    int *ipiv, idx, idx1;
    int i, j, nprobe; 
    int ni[MAX_BLOCKS], ntot, ndim;
    int N, N1, N2;
    double tem;
    

    N = ct.num_blocks;
    for (i = 0; i < N; i++)
    {
        ni[i] = ct.block_dim[i];
    }


    ndim = 0;
    for (i = 0; i < N; i++)
    {
        ndim += ni[i];
    }

    if (ndim != nC)
    {
        printf ("\n %d %d ndim  not equaol to nC in Sgreen_c_noneq_p.c", ndim, nC);
        exit (0);
    }


    ntot = pmo.ntot_low;
    /* allocate matrix and initialization  */
    my_malloc_init( H_tri, ntot, complex double );
 
    matrix_kpoint_center(H_tri, Stri, Htri, creal(ene), ct.kp[pct.kstart].kpt[1], ct.kp[pct.kstart].kpt[2]);

    
    for (nprobe = 0; nprobe < cei.num_probe; nprobe++)
	{	
        N1 = cei.probe_in_block[nprobe];
        N2 = sigma_idx[nprobe];


        for (i = 0; i < pmo.mxllda_cond[N1] * pmo.mxlocc_cond[N1]; i++)
        {
            H_tri[pmo.diag_begin[N1] + i] -= sigma[N2 + i];
        }
    }



  matrix_inverse_anyprobe (H_tri, N, ni, iprobe, Green_C_row, Green_C_col); 
//#if GPU_ENABLED
 // matrix_inverse_anyprobe_cuda (H_tri, N, ni, iprobe, Green_C_row, Green_C_col); 
//#else
//   matrix_inverse_anyprobe_p (H_tri, N, ni, iprobe, Green_C_row, Green_C_col); 
//#endif




    my_free( H_tri );

}
