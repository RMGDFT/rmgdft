/************************** SVN Revision Information **************************
 **    $Id: Sgreen_cond_p.c 985 2008-07-28 19:26:34Z ksaha $    **
******************************************************************************/
 
/*
 *  Calculate the right-up block for transmission calculations.
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "md.h"
#include "pmo.h"


void Sgreen_cond_p (REAL *Htri, REAL *Stri, doublecomplex *sigma_all, int *sigma_idx,
                    REAL eneR, REAL eneI, doublecomplex *green_C, int nC, int iprobe1, int iprobe2)
{
/*   H00, S00: nC * nC real matrix
 *   sigma:  nC * nC complex matrix
 *   green_C  nC * nC complex matrix , output green function
 *   nC: number of states in conductor region
 */

    doublecomplex *H_tri, *green_all;
    /*doublecomplex *H_whole, *H_inv;*/

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
    my_malloc_init( H_tri, ntot, doublecomplex );
    my_malloc_init( ipiv, nC, int );

    /* Construct H = ES - H */
    for (i = 0; i < ntot; i++)
    {
        H_tri[i].r = eneR * Stri[i] - Htri[i] * Ha_eV;
        H_tri[i].i = eneI * Stri[i];
    }


/* --------------------------------------------------- */
	
    /* put the sigma for a probe in the corresponding block 
	   of the Green's matrices  */

    for (nprobe = 0; nprobe < cei.num_probe; nprobe++)
    {
        N1 = cei.probe_in_block[nprobe];
        N2 = sigma_idx[nprobe];
        for (i = 0; i < pmo.mxllda_cond[N1] * pmo.mxlocc_cond[N1]; i++)
        {
            H_tri[pmo.diag_begin[N1] + i].r -= sigma_all[N2 + i].r;
            H_tri[pmo.diag_begin[N1] + i].i -= sigma_all[N2 + i].i;
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

    my_malloc_init( green_all, maxcol * totrow, doublecomplex );

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
            green_C[idx].r = green_all[n_begin1[n1] + idx].r;       
            green_C[idx].i = green_all[n_begin1[n1] + idx].i;       

        }
    }


    time2 = my_crtc ();
    md_timings (matrix_inverse_cond_TIME, (time2 - time1));


    my_free( H_tri );
    my_free(ipiv);


}
