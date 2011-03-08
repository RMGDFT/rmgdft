/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "md.h"
#include "pmo.h"


void sigma_all_energy_point (doublecomplex * sigma_all)
{
    int iprobe, jprobe, idx_delta, j;
    int iene;
    int st1, n2;
    REAL eneR;
    REAL eneI;
    doublecomplex *sigma, *tot, *tott, *g;          

    int idx_sigma, idx_C;
    double time1, time2;
    int nmax, maxrow, maxcol, maxrow2, maxcol2;
    int max_sigma_col, max_sigma_row;
    int t1, t2;

    time1 = my_crtc ();


    nmax = 0;
    maxrow =0;
    maxcol =0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        nmax = max (nmax, lcr[iprobe].num_states);
        maxrow = max( maxrow, pmo.mxllda_lead[iprobe-1]);
        maxcol = max( maxcol, pmo.mxlocc_lead[iprobe-1]);
    }


    my_malloc_init( tot, maxrow * maxcol, doublecomplex);
    my_malloc_init( tott, maxrow * maxcol, doublecomplex);
    my_malloc_init( g, maxrow * maxcol, doublecomplex);

    max_sigma_col = 0;
    max_sigma_row = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        max_sigma_row = pmo.mxllda_cond[idx_C];
        max_sigma_col = pmo.mxlocc_cond[idx_C];
        my_malloc_init( sigma, max_sigma_row * max_sigma_col, doublecomplex);
    }

    /************************/
    int idx, i;



    /*  Calculating the equilibrium term eq. 32 of PRB 65, 165401  */

    idx_sigma = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {

        /*   parallel for the processor on energy grid */
        for (iene = pmo.myblacs; iene < lcr[iprobe].nenergy; iene += pmo.npe_energy)
        {

            eneR = lcr[iprobe].eneR[iene];
            eneI = lcr[iprobe].eneI[iene];


            /* sigma is a complex matrix with dimension ct.num_states * ct.num_states 
             * it sums over all probes
             * tot, tott,  is also a complex matrix, 
             * their memory should be the maximum of probe dimensions, lcr[1,...].num_states
             */


            for (jprobe = 1; jprobe <= cei.num_probe; jprobe++)
            {

                if (eneI > 0.5 )
                {
                    Sgreen_semi_infinite_p (g, lcr[jprobe].H00, lcr[jprobe].H01,
                            lcr[jprobe].S00, lcr[jprobe].S01, eneR, eneI,
                            jprobe);
                }
                else
                {

                    Stransfer_p (tot, tott, lcr[jprobe].H00, lcr[jprobe].H01,
                            lcr[jprobe].S00, lcr[jprobe].S01, eneR, eneI,
                            jprobe);
                    Sgreen_p (tot, tott, lcr[jprobe].H00, lcr[jprobe].H01, lcr[jprobe].S00,
                            lcr[jprobe].S01, eneR, eneI, g, jprobe);

                }

                Sigma_p (sigma, lcr[jprobe].HCL, lcr[jprobe].SCL, eneR, eneI, g, jprobe);

                /*-------------------------------------------------------------------*/

                idx_C = cei.probe_in_block[jprobe - 1];  /* block index */
                for (st1 = 0; st1 < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; st1++)
                {
                    sigma_all[idx_sigma + st1].r = sigma[st1].r;
                    sigma_all[idx_sigma + st1].i = sigma[st1].i;
                }
                idx_sigma += pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C];


            }                   /* end for jprobe */  

        }                       /* end for energy points */

    }                           /* end for iprobe */


    /*  Calculating the non-equilibrium term eq. 33 of PRB 65, 165401  */

    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        j = 0;
        for (idx_delta = 1; idx_delta <= cei.num_probe; idx_delta++)
        {
            if(idx_delta != iprobe)
            {

                /*  parallel for the processor on energy grid */
                for (iene = pmo.myblacs; iene < lcr[iprobe].lcr_ne[j].nenergy_ne; iene += pmo.npe_energy)
                {

                    eneR = lcr[iprobe].lcr_ne[j].eneR_ne[iene];
                    eneI = lcr[iprobe].lcr_ne[j].eneI_ne[iene];


                    /* sigma is a complex matrix with dimension ct.num_states * ct.num_states 
                     * it sums over all probes
                     * tot, tott,  is also a complex matrix, 
                     * their memory should be the maximum of probe dimensions, lcr[1,...].num_states
                     */


                    for (jprobe = 1; jprobe <= cei.num_probe; jprobe++)
                    {

                        if (eneI > 0.5 )
                        {
                            Sgreen_semi_infinite_p (g, lcr[jprobe].H00, lcr[jprobe].H01,
                                    lcr[jprobe].S00, lcr[jprobe].S01, eneR, eneI,
                                    jprobe);
                        }
                        else
                        {

                            Stransfer_p (tot, tott, lcr[jprobe].H00, lcr[jprobe].H01,
                                    lcr[jprobe].S00, lcr[jprobe].S01, eneR, eneI, jprobe);

                            Sgreen_p (tot, tott, lcr[jprobe].H00, lcr[jprobe].H01, lcr[jprobe].S00,
                                    lcr[jprobe].S01, eneR, eneI, g, jprobe);

                        }

                        Sigma_p (sigma, lcr[jprobe].HCL, lcr[jprobe].SCL, eneR, eneI, g, jprobe);


                        idx_C = cei.probe_in_block[jprobe - 1];  /* block index */
                        for (st1 = 0; st1 < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; st1++)
                        {
                            sigma_all[idx_sigma + st1].r = sigma[st1].r;
                            sigma_all[idx_sigma + st1].i = sigma[st1].i;
                        }
                        idx_sigma += pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C];


                    }                         /* jprobe loop ends here */

                }                      /* energy loop ends here */

                j++;
            }            /* if statement ends here */

        }               /* idx_delta loop ends here */

    }        /* iprobe loop ends here */




    my_free(tot);
    my_free(tott);
    my_free(g);
    my_free(sigma);

    time2 = my_crtc ();
    md_timings (SIGMA_ALL_TIME, (time2 - time1));


}
