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


void Sgreen_c_p (rmg_double_t * Htri, rmg_double_t * Stri, complex double * sigma, int * sigma_idx,
                       complex double ene, complex double * Green_C)
{
/*   H00, S00: nC * nC real matrix
 *   sigma:  nC * nC complex matrix
 *   Green_C  nC * nC complex matrix , output green function
 *   nC: number of states in conductor region
 */

    complex double *H_tri;


    int info;
    int i, j, nprobe;
    int ntot, N1, N2; 
    int idx, idx2, ioff, joff;
    

    ntot = pmo.ntot_low;


//  H_tri store the matrix of diagonal, upper offdiag for all blocks
//  then the lower offdiag blocks.

    /* allocate matrix and initialization  */
    my_malloc_init( H_tri, ntot, complex double );
 
    matrix_kpoint_center(H_tri, Stri, Htri, ene, ct.kp[pct.kstart].kpt[1], ct.kp[pct.kstart].kpt[2]);


    /* Construct H = ES - H */
//    for (i = 0; i < ntot; i++)
 //   {
  //      H_tri[i] = ene * Stri[i] - Htri[i] * Ha_eV;
   // }

	
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


//#if GPU_ENABLED
//    matrix_inverse_cuda (H_tri, Green_C);
//#else
    matrix_inverse_p (H_tri, Green_C);
//#endif




    my_free( H_tri );


}

