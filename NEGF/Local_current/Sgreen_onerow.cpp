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

#include "negf_prototypes.h"
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"


void Sgreen_onerow (std::complex<double> *Htri, std::complex<double> * sigma,
                     int *sigma_idx, std::complex<double> * Green_C, int nC,
                     int iprobe)
{
/*   H00, S00: nC * nC real matrix
 *   sigma:  nC * nC complex matrix
 *   Green_C  nC * nC complex matrix , output green function
 *   nC: number of states in conductor region
 */

    std::complex<double> *H_tri;

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
        error_handler ("\n %d %d ndim  not equaol to nC ", ndim, nC);
    }


    ntot = pmo.ntot;
    /* allocate matrix and initialization  */
    my_malloc_init( H_tri, ntot, std::complex<double> );

    
 	/* Construct: H = ES - H */
    for (i = 0; i < ntot; i++)
    {
        H_tri[i] = Htri[i];
    }
    /* put the sigma for a probe in the corresponding block
	   of the Green's matrices  */

    for (nprobe = 0; nprobe < cei.num_probe; nprobe++)
	{	
        N1 = cei.probe_in_block[nprobe];
        N2 = sigma_idx[nprobe];


        for (i = 0; i < pmo.mxllda_cond[N1] * pmo.mxlocc_cond[N1]; i++)
        {
            H_tri[pmo.diag_begin[N1] + i] -= sigma[N2 + i];
        }
    }




   matrix_inverse_anyprobe_p (H_tri, N, ni, iprobe, Green_C); 




    my_free( H_tri );

}
