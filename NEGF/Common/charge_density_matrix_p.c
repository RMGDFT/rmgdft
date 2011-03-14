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

#define 	LDEBUG 	0

void charge_density_matrix_p (doublecomplex * sigma_all)
{
    int iprobe, jprobe, nprobe;
    int iene;
    int st1, idx_sigma, idx_C;
    REAL eneR;
    REAL eneI;
    REAL weightR;
    REAL weightI;
    doublecomplex *green_C, *rho_mn;
    doublecomplex *sigma;
	REAL denominator, numerator, dum, sum, *wmn;
    int nC, nL, i, ntot, *sigma_idx, idx_delta, j;
    doublecomplex *green_C_non;
    int maxrow, maxcol;
    double time1, time2, time3, time4, time5, time6;

    time1 = my_crtc ();

    nL = lcr[1].num_states;
    if (nL != ct.block_dim[0])
    {
        printf (" lcr[1].num_states & ct.block_dim[0] are unequal \n");
    }

    nL = lcr[2].num_states;

    if (nL != ct.block_dim[ct.num_blocks - 1])
    {
        printf (" lcr[2].num_states & ct.block_dim[%d] are unequal \n", ct.num_blocks - 1);
    }

/*  allocate memory for green_C, grenn_C is tri-diagonal */
    ntot = pmo.ntot;
    my_malloc_init( green_C, ntot, doublecomplex );
    my_malloc_init( sigma_idx, cei.num_probe, int );

/*   Calculating the equilibrium term eq. 32 of PRB 65, 165401  */

    time3 = my_crtc ();

    idx_sigma = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        for (st1 = 0; st1 < ntot; st1++)
        {
            lcr[iprobe].density_matrix_tri[st1] = 0.0;

        }
    }

    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        /*   parallel for the processor on energy grid */
        for (iene = pmo.myblacs; iene < lcr[iprobe].nenergy; iene += pmo.npe_energy)
        {

            eneR = lcr[iprobe].eneR[iene];
            eneI = lcr[iprobe].eneI[iene];
            weightR = lcr[iprobe].weightR[iene];
            weightI = lcr[iprobe].weightI[iene];


            /* sigma is a complex matrix with dimension ct.num_states * ct.num_states 
             * it sums over all probes
             * tot, tott, green_tem is also a complex matrix, 
             * their memory should be the maximum of probe dimensions, lcr[1,...].num_states
             */


            for (jprobe = 1; jprobe <= cei.num_probe; jprobe++)
            {
                sigma_idx[jprobe - 1] = idx_sigma;

                idx_C = cei.probe_in_block[jprobe-1];
                idx_sigma += pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C];
            }


            Sgreen_c_p (lcr[0].Htri, lcr[0].Stri, sigma_all, sigma_idx, eneR, eneI, green_C);

            for (st1 = 0; st1 < ntot; st1++)
                lcr[iprobe].density_matrix_tri[st1] +=
                    weightR * green_C[st1].i + weightI * green_C[st1].r;

        }                       /* end for energy points */

        time5 = my_crtc ();
        comm_sums (lcr[iprobe].density_matrix_tri, &ntot, COMM_EN1);
        time6 = my_crtc ();
        md_timings (MPISUM_EQ_TIME, (time6 - time5));


        if (ct.runflag == 111 | ct.runflag == 112 |ct.runflag ==1121)
            break;

    }

    my_free( green_C );

    time4 = my_crtc ();
    md_timings (EQ_PART_TIME, (time4 - time3));

    /* ======================= Non-equilibrium part ===================== */

    if (ct.runflag == 113)
    {

        maxrow = 0;
        maxcol = 0;
        for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
        {
            idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
            maxrow = max(maxrow, pmo.mxllda_cond[idx_C]);
            maxcol = max(maxcol, pmo.mxlocc_cond[idx_C]);
        }

        my_malloc_init( sigma, maxrow * maxcol, doublecomplex ); 


        maxrow = 0;
        maxcol = 0;
        int totrow = 0;
        for(i = 0; i < ct.num_blocks; i++)
        {
            maxrow = max(maxrow, pmo.mxllda_cond[i]);
            maxcol = max(maxcol, pmo.mxlocc_cond[i]);
            totrow += pmo.mxllda_cond[i];
        }

        my_malloc_init( green_C_non, maxcol * totrow, doublecomplex ); 
        my_malloc_init( rho_mn, ntot, doublecomplex );

        /*   Calculating the non-equilibrium term eq. 33 of PRB 65, 165401  */
        for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
        {
            j = 0;	
            for (idx_delta = 1; idx_delta <= cei.num_probe; idx_delta++)
            {
                if(idx_delta != iprobe)
                {

                    for (st1 = 0; st1 < ntot; st1++)
                        lcr[iprobe].lcr_ne[j].density_matrix_ne_tri[st1] = 0.0;

                    /*  parallel for the processor on energy grid */
                    for (iene = pmo.myblacs; iene < lcr[iprobe].lcr_ne[j].nenergy_ne; iene += pmo.npe_energy)
                    {

                        eneR = lcr[iprobe].lcr_ne[j].eneR_ne[iene];
                        eneI = lcr[iprobe].lcr_ne[j].eneI_ne[iene];
                        weightR = lcr[iprobe].lcr_ne[j].weightR_ne[iene];
                        weightI = lcr[iprobe].lcr_ne[j].weightI_ne[iene];

                        for (jprobe = 1; jprobe <= cei.num_probe; jprobe++)
                        {
                            sigma_idx[jprobe - 1] = idx_sigma;
                            idx_C = cei.probe_in_block[jprobe - 1];  /* block index */  
                            idx_sigma += pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C];
                        }


                        nC = ct.num_states;

                        idx_C = cei.probe_in_block[idx_delta - 1];
                        for (st1 = 0; st1 < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; st1++)
                        {
                            sigma[st1].r = sigma_all[sigma_idx[idx_delta - 1] + st1].r;
                            sigma[st1].i = sigma_all[sigma_idx[idx_delta - 1] + st1].i;
                        }

                        Sgreen_c_noneq_p (lcr[0].Htri, lcr[0].Stri, sigma_all, sigma_idx, eneR, eneI,
                                green_C_non, nC, idx_delta);



                        /* sigma is the rhomn in output  */
                        time5 = my_crtc ();

                        rho_munu_p (rho_mn, green_C_non, sigma, idx_delta); 

                        time6 = my_crtc ();
                        md_timings (RHO_MUNU_TIME, (time6 - time5));


                        for (st1 = 0; st1 < ntot; st1++)
                        {

                            lcr[iprobe].lcr_ne[j].density_matrix_ne_tri[st1] +=
                                weightR * rho_mn[st1].r - weightI * rho_mn[st1].i;   /* check -ve sign ??? */

                        }


                    }                   /* end for energy points */

                    comm_sums (lcr[iprobe].lcr_ne[j].density_matrix_ne_tri, &ntot, COMM_EN1);

                    j++;	
                }     /* if statement ends here */	
            }     /* idx_delta loop ends here */
        }      /* iprobe loop ends here */

        my_free( green_C_non );
        my_free( rho_mn );
        my_free( sigma_idx );

        my_free( sigma );

        time3 = my_crtc ();
        md_timings (NONEQ_PART_TIME, (time3 - time4));

        /* ========== Calculation of the density matrix ============= */		

        my_malloc_init( wmn, cei.num_probe, REAL );

        for (st1 = 0; st1 < ntot; st1++)
            lcr[0].density_matrix_tri[st1] = 0.0;

        for (st1 = 0; st1 < ntot; st1++)
        {

            /* Calculates the denominator part for the weight*/
            denominator = 0.0;
            for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
            {
                for (idx_delta = 1; idx_delta < cei.num_probe; idx_delta++)
                {
                    dum = lcr[iprobe].lcr_ne[idx_delta - 1].density_matrix_ne_tri[st1];
                    dum = dum * dum;
                    denominator += dum;

                }  /* idx_delta loop ends here */
            }   /* iprobe loop ends here    */

            denominator = denominator * (cei.num_probe - 1.0); 


            /* Calculates the numerator part for the weight */
            for (nprobe = 1; nprobe <= cei.num_probe; nprobe++)
            {

                numerator = 0.0;
                for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
                {
                    if(iprobe != nprobe)
                    {
                        for (idx_delta = 1; idx_delta < cei.num_probe; idx_delta++)
                        {
                            dum = lcr[iprobe].lcr_ne[idx_delta - 1].density_matrix_ne_tri[st1];
                            dum = dum * dum;
                            numerator += dum;

                        }  /* idx_delta loop ends here */
                    }	 /* if statement ends here  */
                }   /* iprobe loop ends here    */

                wmn[nprobe - 1] = numerator / denominator;
                /* if(nprobe != 1) wmn[nprobe - 1] = 0.0; */
                if(denominator < 0.0000001) wmn[nprobe - 1] = 1.0 / cei.num_probe;  /* Check it */
            }   /* nprobe loop ends here    */


            /*
               if(pct.gridpe == 0) 
               printf (" \n omega %d %f %f %f %f %f \n", st1, wmn[0], wmn[1], wmn[2], wmn[3], wmn[0]+wmn[1]+wmn[2]+wmn[3]);
               printf (" \n omega %d %f %f %f \n", st1, wmn[0], wmn[1], wmn[0]+wmn[1]); 
             */

            /* Finally, calculates density matrix */ 

            for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
            {
                sum = lcr[iprobe].density_matrix_tri[st1];	

                for (idx_delta = 1; idx_delta < cei.num_probe; idx_delta++)
                {

                    sum -= lcr[iprobe].lcr_ne[idx_delta - 1].density_matrix_ne_tri[st1];

                }  
                lcr[0].density_matrix_tri[st1] += wmn[iprobe - 1] * sum;
            }   

            /*		
                    lcr[0].density_matrix_tri[st1] = 
                    wmn[1] * (lcr[1].density_matrix_tri[st1] - lcr[2].lcr_ne[0].density_matrix_ne_tri[st1]) 
                    + wmn[0] * (lcr[2].density_matrix_tri[st1] - lcr[1].lcr_ne[0].density_matrix_ne_tri[st1]);
             */



        } /* st1 loop ends here */
        my_free( wmn );
    }     /* if (ct.runflag == 113) statement ends here */ 


    /* =================== Non-equilibrium part ends here ================= */

    for (st1 = 0; st1 < ntot; st1++)
        lcr[0].density_matrix_tri[st1] /= PI;


    if (cei.num_probe > 4)
        error_handler ("probe > 4");

    time2 = my_crtc ();
    md_timings (CHARGE_DEN_MAT_TIME, (time2 - time1));

}
