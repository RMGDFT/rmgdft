/************************** SVN Revision Information **************************
 **    $Id: Sgreen_c_noneq_p.c 1242 2011-02-02 18:55:23Z luw $    **
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


void Sgreen_c_noneq_p (double *H00, double *S00, doublecomplex * sigma,
                     int *sigma_idx, REAL eneR, REAL eneI, doublecomplex * Green_C, int nC,
                     int iprobe)
{
/*   H00, S00: nC * nC real matrix
 *   sigma:  nC * nC complex matrix
 *   Green_C  nC * nC complex matrix , output green function
 *   nC: number of states in conductor region
 */

    doublecomplex *H_tri, *H_whole, *H_inv;

    int *ipiv, idx, idx1;
    int i, j, nprobe; 
    REAL time1, time2;
    int ni[MAX_BLOCKS], ntot, ndim;
    int N, N1, N2;
    REAL tem;


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
        printf ("\n %d %d ndim  not equaol to nC in Sgreen_c.c", ndim, nC);
        exit (0);
    }


    ntot = pmo.ntot;
    /* allocate matrix and initialization  */
    my_malloc_init( H_tri, ntot, doublecomplex );

    
 	/* Construct: H = ES - H */
    for (i = 0; i < ntot; i++)
    {
        H_tri[i].r = eneR * S00[i] - H00[i] * Ha_eV;
        H_tri[i].i = eneI * S00[i];
    }
    /* put the sigma for a probe in the corresponding block
	   of the Green's matrices  */

    for (nprobe = 0; nprobe < cei.num_probe; nprobe++)
	{	
        N1 = cei.probe_in_block[nprobe];
        N2 = sigma_idx[nprobe];


        for (i = 0; i < pmo.mxllda_cond[N1] * pmo.mxlocc_cond[N1]; i++)
        {
            H_tri[pmo.diag_begin[N1] + i].r -= sigma[N2 + i].r;
            H_tri[pmo.diag_begin[N1] + i].i -= sigma[N2 + i].i;
        }
    }


    time1 = my_crtc ();


   matrix_inverse_anyprobe_p (H_tri, N, ni, iprobe, Green_C); 


    time2 = my_crtc ();
    md_timings (matrix_inverse_lr_TIME, (time2 - time1));


    my_free( H_tri );

}
