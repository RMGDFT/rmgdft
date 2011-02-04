/************************** SVN Revision Information **************************
 **    $Id: Sgreen_c_wang.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "md.h"
#include "pmo.h"


void Sgreen_c_wang (REAL *Htri, REAL *Stri, doublecomplex *sigma_all, int *sigma_idx, 
                    REAL eneR, REAL eneI, doublecomplex *Green_C, int nC)
{
/*   H00, S00: nC * nC real matrix
 *   sigma:  nC * nC complex matrix
 *   Green_C  nC * nC complex matrix , output green function
 *   nC: number of states in conductor region
 */

    doublecomplex *H_tri;


    int info;
    int i, j;
    REAL time1, time2, time3;
    int ni[MAX_BLOCKS], ntot, ndim;
    int N, nprobe, N1, N2;
    REAL tem;

    time3 = my_crtc ();

    N = ct.num_blocks;
    for (i = 0; i < N; i++)
    {
        ni[i] = ct.block_dim[i];
    }

    ntot = pmo.diag_begin[ct.num_blocks-1] + pmo.mxllda_cond[ct.num_blocks-1] * pmo.mxlocc_cond[ct.num_blocks-1];

/*
    ntot = 0;
    ndim = 0;
    for (i = 0; i < N; i++)
    {
        ntot += ni[i] * ni[i];
        ndim += ni[i];
    }

    if (ndim != nC)
    {
        printf ("\n %d %d ndim  not equaol to nC in Sgreen_c.c", ndim, nC);
        exit (0);
    }


    for (i = 1; i < N; i++)
    {
        ntot += ni[i - 1] * ni[i];
    }
*/


    /* allocate matrix and initialization  */
    my_malloc_init( H_tri, ntot, doublecomplex );


    /* Construct H = ES - H */
    for (i = 0; i < ntot; i++)
    {
        H_tri[i].r = eneR * Stri[i] - Htri[i] * Ha_eV;
        H_tri[i].i = eneI * Stri[i];
    }


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

/*
    matrix_inverse_luw (H_tri, Green_C, N, ni); 
*/

    doublecomplex *H_whole, *H_inv;
    int *ipiv;
    my_malloc_init( H_whole, nC * nC, doublecomplex );
    my_malloc_init( H_inv, nC * nC, doublecomplex );
    my_malloc_init( ipiv, nC, int );

    tri_to_whole_complex_p (H_tri, H_whole, N, ni);

    get_inverse_block (H_whole, H_inv, ipiv, nC);

    whole_to_tri_complex (Green_C, H_inv, N, ni);

    my_free( H_whole );
    my_free( H_inv );
    my_free( ipiv );


    time2 = my_crtc ();
    md_timings (matrix_inverse_luw_TIME, (time2 - time1));
    md_timings (GREEN_EQ_TIME, (time2 - time1));


    my_free( H_tri );

}
