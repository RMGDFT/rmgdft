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

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"
#include "GpuAlloc.h"


void Sgreen_c_noneq_p (double *Htri, double *Stri, std::complex<double> * sigma,
                     int *sigma_idx, std::complex<double> ene, std::complex<double> *Green_C, std::complex<double> *Green_C_row, 
                     std::complex<double> *Green_C_col, int nC, int iprobe)
{
/*   H00, S00: nC * nC real matrix
 *   sigma:  nC * nC complex matrix
 *   Green_C  nC * nC complex matrix , output green function
 *   nC: number of states in conductor region
 */

    std::complex<double> *H_tri;

    int i, nprobe; 
    int ni[MAX_BLOCKS], ntot, ndim;
    int N, N1, N2;

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
    H_tri = (std::complex<double> *)RmgMallocHost( ntot *  sizeof(std::complex<double> ));
 
    matrix_kpoint_center(H_tri, Stri, Htri, ene, ct.kp[pct.kstart].kpt[1], ct.kp[pct.kstart].kpt[2]);

    
    for (nprobe = 0; nprobe < cei.num_probe; nprobe++)
	{	
        N1 = cei.probe_in_block[nprobe];
        N2 = sigma_idx[nprobe];


        for (i = 0; i < pmo.mxllda_cond[N1] * pmo.mxlocc_cond[N1]; i++)
        {
            H_tri[pmo.diag_begin[N1] + i] -= sigma[N2 + i];
        }
    }



    matrix_inverse_rowcol (H_tri, iprobe, Green_C, Green_C_row, Green_C_col); 


    RmgFreeHost( H_tri );

}
