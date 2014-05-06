/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <complex.h>

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"


void sigma_all_energy_point (complex double * sigma_all)
{
    int iprobe, jprobe, idx_delta, j;
    int iene;
    int st1, n2;
    complex double *sigma, *tot, *tott, *g;          

    complex double *ch0, *ch01, *ch10;
    double *H10, *S10;
    complex double ene, ctem;

    int idx_sigma, idx_C;
    double time1, time2;
    int  maxrow, maxcol, maxrow2, maxcol2;
    int max_sigma_col, max_sigma_row;
    int t1, t2;
    double one = 1.0, zero = 0.0;
    int ione = 1;
    int *desca, *descb, numst, numstC;

    time1 = my_crtc ();


    maxrow =0;
    maxcol =0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        maxrow = max( maxrow, pmo.mxllda_cond[idx_C]);
        maxcol = max( maxcol, pmo.mxlocc_lead[iprobe-1]);
    }


    my_malloc_init( tot, maxrow * maxcol, complex double);
    my_malloc_init( tott, maxrow * maxcol, complex double);
    my_malloc_init( g, maxrow * maxcol, complex double);
    my_malloc_init( ch0, maxrow * maxcol, complex double);
    my_malloc_init( ch01, maxrow * maxcol, complex double);
    my_malloc_init( ch10, maxrow * maxcol, complex double);
    my_malloc_init( S10, maxrow * maxcol,  double);
    my_malloc_init( H10, maxrow * maxcol,  double);

    max_sigma_col = 0;
    max_sigma_row = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        max_sigma_row = pmo.mxllda_cond[idx_C];
        max_sigma_col = pmo.mxlocc_cond[idx_C];
        my_malloc_init( sigma, max_sigma_row * max_sigma_col, complex double);
    }

    /************************/
    int idx, i;



    /*  Calculating the equilibrium term eq. 32 of PRB 65, 165401  */

    idx_sigma = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {

        if(cei.probe_noneq > 0) iprobe = cei.probe_noneq;
        /*   parallel for the processor on energy grid */
        for (iene = pmo.myblacs; iene < lcr[iprobe].nenergy; iene += pmo.npe_energy)
        {

            ene = lcr[iprobe].ene[iene];


            /* sigma is a complex matrix with dimension ct.num_states * ct.num_states 
             * it sums over all probes
             * tot, tott,  is also a complex matrix, 
             * their memory should be the maximum of probe dimensions, lcr[1,...].num_states
             */


            for (jprobe = 1; jprobe <= cei.num_probe; jprobe++)
            {

                desca = &pmo.desc_lead[ (jprobe-1) * DLEN];

                numst = lcr[jprobe].num_states;

                PDTRAN(&numst, &numst, &one, lcr[jprobe].S01, &ione, &ione, desca,
                        &zero, S10, &ione, &ione, desca);
                PDTRAN(&numst, &numst, &one, lcr[jprobe].H01, &ione, &ione, desca,
                        &zero, H10, &ione, &ione, desca);


                idx = pmo.mxllda_lead[jprobe-1] * pmo.mxlocc_lead[jprobe-1];
                for (i = 0; i < idx; i++)
                {
                    ch0[i] = ene * lcr[jprobe].S00[i] - Ha_eV * lcr[jprobe].H00[i];
                    ch01[i] = ene * lcr[jprobe].S01[i] - Ha_eV * lcr[jprobe].H01[i];
                    ch10[i] = ene * S10[i] - Ha_eV * H10[i];
                }

//#if GPU_ENABLED
//                    Sgreen_cuda (g, ch0, ch01, ch10, jprobe);

//#else

                if (cimag(ene) >0.5 )
                {
                    Sgreen_semi_infinite_p (g, ch0, ch01, jprobe);
                }
                else
                {

                    Stransfer_p (tot, tott, ch0, ch01, ch10, jprobe);
                    Sgreen_p (tot, tott, ch0, ch01, g, jprobe);

                }
//#endif

                
                idx_C = cei.probe_in_block[jprobe - 1];  /* block index */
                idx = pmo.mxllda_cond[idx_C] * pmo.mxlocc_lead[jprobe-1];
                for (i = 0; i < idx; i++)
                {

                    ch01[i] = ene * lcr[jprobe].SCL[i] - Ha_eV * lcr[jprobe].HCL[i];
                }
                desca = &pmo.desc_cond_lead[ (idx_C + (jprobe-1) * ct.num_blocks) * DLEN];
                descb = &pmo.desc_lead_cond[ (idx_C + (jprobe-1) * ct.num_blocks) * DLEN];
                numst = lcr[jprobe].num_states;
                numstC = ct.block_dim[idx_C];


                PDTRAN(&numst, &numstC, &one, lcr[jprobe].SCL, &ione, &ione, desca,
                        &zero, S10, &ione, &ione, descb);
                PDTRAN(&numst, &numstC, &one, lcr[jprobe].HCL, &ione, &ione, desca,
                        &zero, H10, &ione, &ione, descb);
                idx = pmo.mxllda_lead[jprobe -1] * pmo.mxlocc_cond[idx_C];
                for (i = 0; i < idx; i++)
                {
                    ch10[i] = ene * S10[i] - Ha_eV * H10[i];
                }

                Sigma_p (sigma, ch0, ch01, ch10, g, jprobe);

                /*-------------------------------------------------------------------*/

                idx_C = cei.probe_in_block[jprobe - 1];  /* block index */
                for (st1 = 0; st1 < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; st1++)
                {
                    sigma_all[idx_sigma + st1] = sigma[st1];
                }
                idx_sigma += pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C];


            }                   /* end for jprobe */  

        }                       /* end for energy points */

        if(cei.probe_noneq > 0) break;

    }                           /* end for iprobe */


    /*  Calculating the non-equilibrium term eq. 33 of PRB 65, 165401  */

    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {

        if(cei.probe_noneq >0) iprobe = cei.probe_noneq;
        j = 0;
        for (idx_delta = 1; idx_delta <= cei.num_probe; idx_delta++)
        {
            if(idx_delta != iprobe)
            {

                /*  parallel for the processor on energy grid */
                for (iene = pmo.myblacs; iene < lcr[iprobe].lcr_ne[j].nenergy_ne; iene += pmo.npe_energy)
                {

                    ene = lcr[iprobe].lcr_ne[j].ene_ne[iene];


                    /* sigma is a complex matrix with dimension ct.num_states * ct.num_states 
                     * it sums over all probes
                     * tot, tott,  is also a complex matrix, 
                     * their memory should be the maximum of probe dimensions, lcr[1,...].num_states
                     */


                    for (jprobe = 1; jprobe <= cei.num_probe; jprobe++)
                    {

                        desca = &pmo.desc_lead[ (jprobe-1) * DLEN];

                        numst = lcr[jprobe].num_states;

                        PDTRAN(&numst, &numst, &one, lcr[jprobe].S01, &ione, &ione, desca,
                                &zero, S10, &ione, &ione, desca);
                        PDTRAN(&numst, &numst, &one, lcr[jprobe].H01, &ione, &ione, desca,
                                &zero, H10, &ione, &ione, desca);


                        idx = pmo.mxllda_lead[jprobe-1] * pmo.mxlocc_lead[jprobe-1];
                        for (i = 0; i < idx; i++)
                        {
                            ch0[i] = ene * lcr[jprobe].S00[i] - Ha_eV * lcr[jprobe].H00[i];
                            ch01[i] = ene * lcr[jprobe].S01[i] - Ha_eV * lcr[jprobe].H01[i];
                            ch10[i] = ene * S10[i] - Ha_eV * H10[i];
                        }

//#if GPU_ENABLED
//                    Sgreen_cuda (g, ch0, ch01, ch10, jprobe);

//#else
                        Stransfer_p (tot, tott, ch0, ch01, ch10, jprobe);
                        Sgreen_p (tot, tott, ch0, ch01, g, jprobe);
//#endif



                        idx_C = cei.probe_in_block[jprobe - 1];  /* block index */
                        idx = pmo.mxllda_cond[idx_C] * pmo.mxlocc_lead[jprobe-1];
                        for (i = 0; i < idx; i++)
                        {
                            ch01[i] = ene * lcr[jprobe].SCL[i] - Ha_eV * lcr[jprobe].HCL[i];
                        }
                        desca = &pmo.desc_cond_lead[ (idx_C + (jprobe-1) * ct.num_blocks) * DLEN];
                        descb = &pmo.desc_lead_cond[ (idx_C + (jprobe-1) * ct.num_blocks) * DLEN];
                        numst = lcr[jprobe].num_states;
                        numstC = ct.block_dim[idx_C];


                        PDTRAN(&numst, &numstC, &one, lcr[jprobe].SCL, &ione, &ione, desca,
                                &zero, S10, &ione, &ione, descb);
                        PDTRAN(&numst, &numstC, &one, lcr[jprobe].HCL, &ione, &ione, desca,
                                &zero, H10, &ione, &ione, descb);
                        idx = pmo.mxllda_lead[jprobe -1] * pmo.mxlocc_cond[idx_C];
                        for (i = 0; i < idx; i++)
                        {
                            ch10[i] = ene * S10[i] - Ha_eV * H10[i];
                        }

                        Sigma_p (sigma, ch0, ch01, ch10, g, jprobe);


                        idx_C = cei.probe_in_block[jprobe - 1];  /* block index */
                        for (st1 = 0; st1 < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; st1++)
                        {
                            sigma_all[idx_sigma + st1] = sigma[st1];
                        }
                        idx_sigma += pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C];


                    }                         /* jprobe loop ends here */

                }                      /* energy loop ends here */

                j++;
            }            /* if statement ends here */

        }               /* idx_delta loop ends here */

        if(cei.probe_noneq >0) break;

    }        /* iprobe loop ends here */




    my_free(tot);
    my_free(tott);
    my_free(g);
    my_free(sigma);
    my_free(ch0);
    my_free(ch01);
    my_free(ch10);

    time2 = my_crtc ();
    rmg_timings (SIGMA_ALL_TIME, (time2 - time1));


}
