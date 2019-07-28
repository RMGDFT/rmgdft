#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
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


// calculating sigma for jprobe at energy ene

void sigma_one_energy_point (std::complex<double> *sigma, int jprobe, std::complex<double> ene, 
double kvecy, double kvecz, std::complex<double> *work)
{
    int iprobe, idx_delta, j;
    int st1, n2;
    std::complex<double> *g;          
    std::complex<double> *ch0, *ch01, *ch10;
    std::complex<double> *S00, *H00, *H10, *S10, *H01, *S01, *HCL, *SCL;

    int idx_sigma, idx_C;
    int  maxrow, maxcol, maxrow2, maxcol2;
    std::complex<double> one = 1.0, zero = 0.0;
    int ione = 1;
    int *desca, *descb, numst, numstC;



    maxrow =0;
    maxcol =0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        maxrow = rmg_max( maxrow, pmo.mxllda_cond[idx_C]);
        maxcol = rmg_max( maxcol, pmo.mxlocc_lead[iprobe-1]);
    }


    n2 = maxrow * maxcol;
    g    = work + 0*n2;
    ch0  = work + 1*n2;
    ch01 = work + 2*n2;
    ch10 = work + 3*n2;
    S00  = work + 4*n2;
    H00  = work + 5*n2;
    S10  = work + 6*n2;
    H10  = work + 7*n2;
    S01  = work + 8*n2;
    H01  = work + 9*n2;
    SCL  = work +10*n2;
    HCL  = work +11*n2;


    int idx, i;


    matrix_kpoint_lead(S00, H00, S01, H01, SCL, HCL,  kvecy, kvecz, jprobe);
    desca = &pmo.desc_lead[ (jprobe-1) * DLEN];

    numst = lcr[jprobe].num_states;

    pztranc(&numst, &numst, &one, S01, &ione, &ione, desca,
            &zero, S10, &ione, &ione, desca);
    pztranc(&numst, &numst, &one, H01, &ione, &ione, desca,
            &zero, H10, &ione, &ione, desca);


    idx = pmo.mxllda_lead[jprobe-1] * pmo.mxlocc_lead[jprobe-1];
 //   for (i = 0; i < idx; i++)
 //   {
 //       ch0[i]  = std::real(ene) * S00[i] - Ha_eV * H00[i];
 //       ch01[i] = std::real(ene) * S01[i] - Ha_eV * H01[i];
 //       ch10[i] = std::real(ene) * S10[i] - Ha_eV * H10[i];
 //   }

    for (i = 0; i < idx; i++)
    {
        ch0[i]  = ene * S00[i] - Ha_eV * H00[i];
        ch01[i] = ene * S01[i] - Ha_eV * H01[i];
        ch10[i] = ene * S10[i] - Ha_eV * H10[i];
    }


    if (std::imag(ene) <0.5 )
    {

        //KrylovSigma_c(numst, ch0, ch10, ch01,sigma, 0.01);
        //return;
        green_lead(ch0, ch01, ch10, g, jprobe);

    }    
    else
    {
        Sgreen_semi_infinite_p (g, ch0, ch01, ch10, jprobe);
    }
    //    else
    //    {
    //
    //
    //    }
    //#endif


    idx_C = cei.probe_in_block[jprobe - 1];  /* block index */
    idx = pmo.mxllda_cond[idx_C] * pmo.mxlocc_lead[jprobe-1];
    for (i = 0; i < idx; i++)
    {

        ch01[i] = ene * SCL[i] - Ha_eV * HCL[i];
    }
    desca = &pmo.desc_cond_lead[ (idx_C + (jprobe-1) * ct.num_blocks) * DLEN];
    descb = &pmo.desc_lead_cond[ (idx_C + (jprobe-1) * ct.num_blocks) * DLEN];
    numst = lcr[jprobe].num_states;
    numstC = ct.block_dim[idx_C];


    pztranc(&numst, &numstC, &one, SCL, &ione, &ione, desca,
            &zero, S10, &ione, &ione, descb);
    pztranc(&numst, &numstC, &one, HCL, &ione, &ione, desca,
            &zero, H10, &ione, &ione, descb);
    idx = pmo.mxllda_lead[jprobe -1] * pmo.mxlocc_cond[idx_C];
    for (i = 0; i < idx; i++)
    {
        ch10[i] = ene * S10[i] - Ha_eV * H10[i];
    }

    Sigma_p (sigma, ch0, ch01, ch10, g, jprobe);

}
