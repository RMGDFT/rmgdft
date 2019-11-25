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


void sigma_all_energy_point (std::complex<double> * sigma_all, double kvecy, double kvecz)
{
    int iprobe, jprobe, idx_delta, j;
    int iene;
    int st1;
    std::complex<double> *sigma;          

    std::complex<double> *work;
    std::complex<double> ene, ctem;

    int idx_sigma, idx_C;
    int  maxrow, maxcol;
    int max_sigma_col, max_sigma_row;



    maxrow =0;
    maxcol =0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        maxrow = rmg_max( maxrow, pmo.mxllda_cond[idx_C]);
        maxcol = rmg_max( maxcol, pmo.mxlocc_lead[iprobe-1]);
    }


    my_malloc_init( work, 12*maxrow * maxcol, std::complex<double>);

    max_sigma_col = 0;
    max_sigma_row = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        max_sigma_row = rmg_max(max_sigma_row, pmo.mxllda_cond[idx_C]);
        max_sigma_col = rmg_max(max_sigma_col, pmo.mxlocc_cond[idx_C]);
    }

    my_malloc_init( sigma, max_sigma_row * max_sigma_col, std::complex<double>);

    /************************/

    /*  Calculating the equilibrium term eq. 32 of PRB 65, 165401  */

    idx_sigma = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {

        if(cei.probe_noneq > 0) iprobe = cei.probe_noneq;
        /*   parallel for the processor on energy grid */
        for (iene = pmo.myblacs; iene < lcr[iprobe].nenergy; iene += pmo.npe_energy)
        {

            ene = lcr[iprobe].ene[iene];



            for (jprobe = 1; jprobe <= cei.num_probe; jprobe++)
            {

                sigma_one_energy_point(sigma, jprobe, ene, kvecy, kvecz, work);

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



                    for (jprobe = 1; jprobe <= cei.num_probe; jprobe++)
                    {

                        sigma_one_energy_point(sigma, jprobe, ene, kvecy, kvecz, work);


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




    my_free(work);
    my_free(sigma);

}
