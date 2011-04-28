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
#include <complex.h>

#include "md.h"
#include "pmo.h"


void Sgreen_cond_p (complex double *H_tri, complex double *sigma_all, int *sigma_idx,
                    complex double *green_C, int nC, int iprobe1, int iprobe2)
{
/*   H00, S00: nC * nC real matrix
 *   sigma:  nC * nC complex matrix
 *   green_C  nC * nC complex matrix , output green function
 *   nC: number of states in conductor region
 */

    complex double *green_all;
    /*complex double *H_whole, *H_inv;*/

    int info;
    int i, nprobe;
    REAL time1, time2;
    int ni[MAX_BLOCKS], ntot ;
    int N, N1, N2;
    int *ipiv;
    int j, idx, idx1;
    int n1, n2, maxcol, totrow, *n_begin1;

    N = ct.num_blocks;
    for (i = 0; i < N; i++)
    {
        ni[i] = ct.block_dim[i];
    }

    ntot = pmo.diag_begin[ct.num_blocks-1] + pmo.mxllda_cond[ct.num_blocks-1] * pmo.mxlocc_cond[ct.num_blocks-1];


    /* allocate matrix and initialization  */
    my_malloc_init( ipiv, nC, int );


/* --------------------------------------------------- */
	
    /* put the sigma for a probe in the corresponding block 
	   of the Green's matrices  */

    for (nprobe = 0; nprobe < cei.num_probe; nprobe++)
    {
        N1 = cei.probe_in_block[nprobe];
        N2 = sigma_idx[nprobe];
        for (i = 0; i < pmo.mxllda_cond[N1] * pmo.mxlocc_cond[N1]; i++)
        {
            H_tri[pmo.diag_begin[N1] + i] -= sigma_all[N2 + i];
        }
    }
    

    time1 = my_crtc ();
/* ===================================================== */

    maxcol = 0;
    totrow = 0;
    for(i = 0; i < ct.num_blocks; i++)
    {
        maxcol = max(maxcol, pmo.mxlocc_cond[i]);
        totrow += pmo.mxllda_cond[i];
    }

    my_malloc_init( green_all, maxcol * totrow, complex double );

    my_malloc( n_begin1, ct.num_blocks, int );
    n_begin1[0] = 0;
    for (i = 1; i < ct.num_blocks; i++)
    {
        n_begin1[i] = n_begin1[i - 1] + pmo.mxllda_cond[i - 1] * maxcol;
    }

    matrix_inverse_anyprobe_p (H_tri, N, ni, iprobe2, green_all);

    n1 = cei.probe_in_block[iprobe1 - 1];
    n2 = cei.probe_in_block[iprobe2 - 1];

    for(i =0; i < pmo.mxlocc_cond[n2]; i++)
    {
        for(j =0; j < pmo.mxllda_cond[n1]; j++)
        {
            idx = j + i * pmo.mxllda_cond[n1];
            green_C[idx] = green_all[n_begin1[n1] + idx];       

        }
    }


    time2 = my_crtc ();
    md_timings (matrix_inverse_cond_TIME, (time2 - time1));


    my_free(ipiv);
    my_free(green_all);


}
