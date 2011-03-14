/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

/*
 *
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "md.h"
#include "pmo.h"

double pmo_trace(doublecomplex*, int*);


void get_cond_frommatrix ()
{
    int iprobe, iprobe1, iprobe2;
    int iene, icond;
    REAL eneR, eneI;
    doublecomplex *tot, *tott, *g;
    doublecomplex *green_C;
    doublecomplex *temp_matrix1, *temp_matrix2;
    doublecomplex *Gamma1, *Gamma2, *sigma, *sigma_all;
    REAL de, emin, emax, E_imag, KT, current;
    REAL *ener1, *cond;
    int ntot, ndim, nC, idx_C, *sigma_idx;
    REAL cons, EF1, EF2, f1, f2;


    doublecomplex alpha, beta;
    int i, j, idx, E_POINTS, kpoint;
    char fcd_n = 'N', fcd_c = 'C', newname[100];
    FILE *file;
    int ione =1, *desca, *descb, *descc, *descd;
    int n1, n2, nC_1, nC_2, nC_11, nC_22, nC_max;
    int idx1, idx2;


/*=============== Reading input and then print them ==============*/ 

    read_cond_input (&emin, &emax, &E_POINTS, &E_imag, &KT, &kpoint);
    de = (emax - emin) / (E_POINTS - 1);


    if (pct.gridpe == 0)
    {
        printf ("\n transmission calculations from known matrix \n");
        for (idx = 0; idx < ct.num_cond_curve; idx++)
        {
            printf ("Calculating transmission from probe %d to %d \n", 
                     ct.cond_probe1[idx], ct.cond_probe2[idx]);
        }	
        printf ("ct.num_states     = %d \n", ct.num_states);
        for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
        {
		 	printf ("lcr[iprobe].num_states = %d \n", lcr[iprobe].num_states);
        }	
        printf ("num of blocks     = %d \n", ct.num_blocks);
        printf ("blocks dim        =   ");
        for (idx = 0; idx < ct.num_blocks; idx++)
            printf (" %d ", ct.block_dim[idx]);
        printf ("\n");

        printf ("enengy from %f to %f with %d points\n", emin, emax, E_POINTS);
        printf ("small imaginary part = %f \n", E_imag);
        printf ("KT = %f eV\n", KT);
    }

/*======================== Reading Matrices ===============================*/

    ntot = 0;
    ndim = 0;
    for (i = 0; i < ct.num_blocks; i++)
    {
        ntot += pmo.mxllda_cond[i] * pmo.mxlocc_cond[i];
        ndim += ct.block_dim[i];
    }
    for (i = 1; i < ct.num_blocks; i++)
    {
        ntot += pmo.mxllda_cond[i-1] * pmo.mxlocc_cond[i];
    }
    if (ndim != ct.num_states)
    {
        printf (" %d %d ndim not equal to nC in get_cond_frommatrix\n", ndim, ct.num_states);
		exit (0);
    }


    my_malloc_init( lcr[0].Htri, ntot, REAL );
    my_malloc_init( lcr[0].Stri, ntot, REAL );

	for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
	{	
		idx = pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];
		my_malloc_init( lcr[iprobe].H00, idx, REAL );
		my_malloc_init( lcr[iprobe].S00, idx, REAL );
		my_malloc_init( lcr[iprobe].H01, idx, REAL );
		my_malloc_init( lcr[iprobe].S01, idx, REAL );
    }

	for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
	{	
		i = cei.probe_in_block[iprobe - 1];
		idx = pmo.mxllda_cond[i] * pmo.mxlocc_lead[iprobe-1];
		my_malloc_init( lcr[iprobe].HCL, idx, REAL );
		my_malloc_init( lcr[iprobe].SCL, idx, REAL );
    }

    read_matrix_pp();

/*=================== Allocate memory for sigma ======================*/

    idx = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        idx = max(idx, pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]);
    }
    my_malloc_init( sigma, idx, doublecomplex );


    my_malloc_init( sigma_idx, cei.num_probe, int ); 

    idx = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        sigma_idx[iprobe - 1] = idx;
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        idx += pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C];
    }

    my_malloc_init( sigma_all, idx, doublecomplex );


/*============== Allocate memory for tot, tott, g ====================*/

    idx = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx = max(idx, pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1]);
    }
 
    my_malloc_init( tot,  idx, doublecomplex );
    my_malloc_init( tott, idx, doublecomplex );
    my_malloc_init( g,    idx, doublecomplex );

/*===================================================================*/

    my_malloc_init( ener1, E_POINTS, REAL );
    my_malloc_init( cond, E_POINTS, REAL );

    alpha.r = 1.0;
    alpha.i = 0.0;
    beta.r = 0.0;
    beta.i = 0.0;

    nC = ct.num_states;

/*===================================================================*/

    for (icond = 0; icond < ct.num_cond_curve; icond++)
    {

        iprobe1 = ct.cond_probe1[icond];
        n1 = cei.probe_in_block[iprobe1 - 1];
        nC_1 = ct.block_dim[n1];
        nC_11 = pmo.mxllda_cond[n1] * pmo.mxlocc_cond[n1];
 
        iprobe2 = ct.cond_probe2[icond];
        n2 = cei.probe_in_block[iprobe2 - 1];
        nC_2 = ct.block_dim[n2];
        nC_22 = pmo.mxllda_cond[n2] * pmo.mxlocc_cond[n2];

        my_malloc_init( Gamma1, nC_11, doublecomplex ); 
        my_malloc_init( Gamma2, nC_22, doublecomplex ); 

        nC_max = max(nC_11, nC_22);

        my_malloc_init( green_C, nC_max, doublecomplex); 
        my_malloc_init( temp_matrix1, nC_max, doublecomplex );
        my_malloc_init( temp_matrix2, nC_max, doublecomplex );

/*===================================================================*/


        for (iene = 0; iene < E_POINTS; iene++)
        {
            ener1[iene] = emin + iene * de;
            cond[iene] = 0.0;
        }


        for (iene = pmo.myblacs; iene < E_POINTS; iene += pmo.npe_energy)
        {

            eneR = emin + iene * de;
            eneI = E_imag;


            for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
            {
                Stransfer_p (tot, tott, lcr[iprobe].H00, lcr[iprobe].H01, 
                        lcr[iprobe].S00, lcr[iprobe].S01, eneR, eneI, iprobe);
    
                Sgreen_p (tot, tott, lcr[iprobe].H00, lcr[iprobe].H01, lcr[iprobe].S00,
                        lcr[iprobe].S01, eneR, eneI, g, iprobe);

                Sigma_p (sigma, lcr[iprobe].HCL, lcr[iprobe].SCL, eneR, eneI, g, iprobe);

                idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
                for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                {
                    sigma_all[sigma_idx[iprobe - 1] + i].r = sigma[i].r;
                    sigma_all[sigma_idx[iprobe - 1] + i].i = sigma[i].i;
                }

                if(iprobe == iprobe1)
                {
                    for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                    {
                        Gamma1[i].r = -2.0 * sigma[i].i;
                        Gamma1[i].i =  0.0;
                    }
                }

                if(iprobe == iprobe2)
                {
                    for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                    {
                        Gamma2[i].r = -2.0 * sigma[i].i;
                        Gamma2[i].i =  0.0;
                    }
                }

            }  /*  end for iprobe */




            Sgreen_cond_p (lcr[0].Htri, lcr[0].Stri, sigma_all, sigma_idx, eneR, eneI, 
                           green_C, nC, iprobe1, iprobe2);



            desca = &pmo.desc_cond[( n1 + n1 * ct.num_blocks) * DLEN];   /* nC_1 * nC_1 matrix */
            descb = &pmo.desc_cond[( n2 + n2 * ct.num_blocks) * DLEN];   /* nC_2 * nC_2 matrix */
            descc = &pmo.desc_cond[( n2 + n1 * ct.num_blocks) * DLEN];   /* nC_2 * nC_1 matrix */
            descd = &pmo.desc_cond[( n1 + n2 * ct.num_blocks) * DLEN];   /* nC_1 * nC_2 matrix */



            /* Gamma(C_1, C_1) * G(C_1, C_2) = temp(C_1, C_2) */
            PZGEMM (&fcd_n, &fcd_n, &nC_1, &nC_2, &nC_1, &alpha, Gamma1, &ione, &ione, desca,
                    green_C, &ione, &ione, descd, &beta, temp_matrix1, &ione, &ione, descd);

            /* temp(C_1, C_2) * Gamma(C_2, C_2) = temp2(C_1, C_2) */
            PZGEMM (&fcd_n, &fcd_n, &nC_1, &nC_2, &nC_2, &alpha, temp_matrix1, &ione, &ione, descd,
                    Gamma2, &ione, &ione, descb, &beta, temp_matrix2, &ione, &ione, descd);

            /* temp2(C_1, C_2) * G(C_2, C_1) = temp(C_1, C_1) */
            PZGEMM (&fcd_n, &fcd_c, &nC_1, &nC_1, &nC_2, &alpha, temp_matrix2, &ione, &ione, descd,
                    green_C, &ione, &ione, descd, &beta, temp_matrix1, &ione, &ione, desca);

            /* desca= &pmo.desc_lead[0]; */
            cond[iene] = pmo_trace(temp_matrix1, desca);

            ener1[iene] = eneR;


            /* printf (" condcond eneR, G= %f %f \n ", eneR, cond[iene]); */


        }                           /*  end for iene */

        my_barrier ();
        iene = E_POINTS;
        global_sums (cond, &iene);

        if (pct.gridpe == 0)
        {
            sprintf(newname, "%s%d%d%s", "cond_", iprobe1, iprobe2, ".dat");
            file = fopen (newname, "w");
            for (iene = 0; iene < E_POINTS; iene++)
                fprintf (file, " %f %22.12f\n", ener1[iene], cond[iene]);
            fclose (file);
        }

        /* calculating the current */

        current = 0.0;
        EF1 = lcr[iprobe1].bias;
        EF2 = lcr[iprobe2].bias;


        for (iene = 1; iene < E_POINTS; iene++)
        {
            f1 = 1.0 / (1.0 + exp ((EF1 - ener1[iene]) / KT));
            f2 = 1.0 / (1.0 + exp ((EF2 - ener1[iene]) / KT));
            current += (f2 - f1) * cond[iene] * (ener1[iene] - ener1[iene - 1]);
        }

        current *= 2.0 * e_C * e_C / h_SI;

        if (pct.gridpe == 0)
        {
            printf ("\n bias = %f eV    current = %e microA\n", EF1 - EF2, current * 1e6);
        }



        my_free(Gamma1);
        my_free(Gamma2);
        my_free(green_C);
        my_free(temp_matrix1);
        my_free(temp_matrix2);

    } /* ends of icond (num_cond_curve) loop */

    my_free(ener1);
    my_free(cond);

    my_free(tot);
    my_free(tott);
    my_free(g);
    my_free(sigma);
    my_free(sigma_all);
    my_free(sigma_idx);

    my_free(lcr[0].Htri);
    my_free(lcr[0].Stri);
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
	{	
		my_free(lcr[iprobe].H00);
		my_free(lcr[iprobe].S00);
		my_free(lcr[iprobe].H01);
		my_free(lcr[iprobe].S01);
		my_free(lcr[iprobe].HCL);
		my_free(lcr[iprobe].SCL);
    }

    my_barrier ();

}
