#include "negf_prototypes.h"
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
#include <complex>

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"


void Sgreen_cond_p (std::complex<double> *H_tri, std::complex<double> *G_tri, std::complex<double> *sigma_all, int *sigma_idx,
                    std::complex<double> *green_C, int nC, int iprobe1, int iprobe2)
{
/*   H00, S00: nC * nC real matrix
 *   sigma:  nC * nC complex matrix
 *   green_C  nC * nC complex matrix , output green function
 *   nC: number of states in conductor region
 */

    int info;
    int i, nprobe;
    int ni[MAX_BLOCKS] ;
    int N, N1, N2;
    int j, idx, idx1;
    int n1, n2, maxcol;

    N = ct.num_blocks;
    for (i = 0; i < N; i++)
    {
        ni[i] = ct.block_dim[i];
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
            H_tri[pmo.diag_begin[N1] + i] -= sigma_all[N2 + i];
        }
    }
    

/* ===================================================== */



    n1 = cei.probe_in_block[iprobe1 - 1];
    n2 = cei.probe_in_block[iprobe2 - 1];
    matrix_inverse_blocknm_Gauss (H_tri, G_tri, n1, n2, green_C);


    /*  Green_C store the (n2,n1) block of Green function */






}
