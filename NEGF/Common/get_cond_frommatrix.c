/*
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
#include <complex.h>

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"

double pmo_trace(complex double*, int*);


void get_cond_frommatrix ()
{
    int iprobe, iprobe1, iprobe2;
    int iene, icond;
    complex double * H_tri,*g;
    complex double *green_C;
    complex double *temp_matrix1, *temp_matrix2;
    complex double *Gamma1, *Gamma2, *sigma;
    double de, emin, emax, E_imag, KT, current;
    double *ener1, *cond;
    int ntot, ndim, nC, idx_C, *sigma_idx;
    double cons, EF1, EF2, f1, f2;

    complex double *ch0, *ch01, *ch10;
    double *H10, *S10;
    complex double ene, ctem;

    complex double alpha, beta;
    double one, zero;
    int i, j, idx, E_POINTS, kpoint[3];
    char fcd_n = 'N', fcd_c = 'C', newname[MAX_PATH];
    FILE *file;
    int ione =1, *desca, *descb, *descc, *descd;
    int n1, n2, nC_1, nC_2, nC_11, nC_22, nC_max;
    int idx1, idx2;
    int numst, numstC;


/*=============== Reading input and then print them ==============*/ 

    read_cond_input (&emin, &emax, &E_POINTS, &E_imag, &KT, kpoint);
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


    my_malloc_init( H_tri, ntot, complex double );
    my_malloc_init( lcr[0].Htri, ntot, double );
    my_malloc_init( lcr[0].Stri, ntot, double );

	for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
	{	
		idx = pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];
		my_malloc_init( lcr[iprobe].H00, idx, double );
		my_malloc_init( lcr[iprobe].S00, idx, double );
		my_malloc_init( lcr[iprobe].H01, idx, double );
		my_malloc_init( lcr[iprobe].S01, idx, double );
    }

	for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
	{	
		i = cei.probe_in_block[iprobe - 1];
		idx = pmo.mxllda_cond[i] * pmo.mxlocc_lead[iprobe-1];
		my_malloc_init( lcr[iprobe].HCL, idx, double );
		my_malloc_init( lcr[iprobe].SCL, idx, double );
    }

    read_matrix_pp();

/*=================== Allocate memory for sigma ======================*/

    idx = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        idx = rmg_max(idx, pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]);
    }
    my_malloc_init( sigma, idx, complex double );


    my_malloc_init( sigma_idx, cei.num_probe, int ); 

    idx = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        sigma_idx[iprobe - 1] = idx;
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        idx += pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C];
    }

    my_malloc_init( sigma_all, idx, complex double );


/*============== Allocate memory for tot, tott, g ====================*/

    idx = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        idx = rmg_max(idx, pmo.mxllda_cond[idx_C] * pmo.mxlocc_lead[iprobe-1]);
    }
 
    my_malloc_init( g,    idx, complex double );
    my_malloc_init( ch0,  idx, complex double );
    my_malloc_init( ch01,  idx, complex double );
    my_malloc_init( ch10,  idx, complex double );

    my_malloc_init( H10,    idx, double );
    my_malloc_init( S10,    idx, double );
/*===================================================================*/

    my_malloc_init( ener1, E_POINTS, double );
    my_malloc_init( cond, E_POINTS, double );

    alpha = 1.0;
    beta = 0.0;
    one = 1.0; 
    zero = 0.0;

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

        my_malloc_init( Gamma1, nC_11, complex double ); 
        my_malloc_init( Gamma2, nC_22, complex double ); 

        nC_max = rmg_max(nC_11, nC_22);

        my_malloc_init( green_C, nC_max, complex double); 
        my_malloc_init( temp_matrix1, nC_max, complex double );
        my_malloc_init( temp_matrix2, nC_max, complex double );

/*===================================================================*/


        for (iene = 0; iene < E_POINTS; iene++)
        {
            ener1[iene] = emin + iene * de;
            cond[iene] = 0.0;
        }


        for (iene = pmo.myblacs; iene < E_POINTS; iene += pmo.npe_energy)
        {

            
            ene = emin + iene * de + I * E_imag;


            for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
            {

                desca = &pmo.desc_lead[ (iprobe-1) * DLEN];

                numst = lcr[iprobe].num_states;

                PDTRAN(&numst, &numst, &one, lcr[iprobe].S01, &ione, &ione, desca, 
                        &zero, S10, &ione, &ione, desca); 
                PDTRAN(&numst, &numst, &one, lcr[iprobe].H01, &ione, &ione, desca, 
                        &zero, H10, &ione, &ione, desca); 

                idx = pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];
                for (i = 0; i < idx; i++)
                {
                    ch0[i] = ene * lcr[iprobe].S00[i] - Ha_eV * lcr[iprobe].H00[i];
                    ch01[i] = ene * lcr[iprobe].S01[i] - Ha_eV * lcr[iprobe].H01[i];
                    ch10[i] = ene * S10[i] - Ha_eV * H10[i];
                }


                green_lead (ch0, ch01, ch10, g, iprobe);

                idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
                idx = pmo.mxllda_cond[idx_C] * pmo.mxlocc_lead[iprobe-1];
                for (i = 0; i < idx; i++)
                {
                    ch01[i] = ene * lcr[iprobe].SCL[i] - Ha_eV * lcr[iprobe].HCL[i];
                }

                desca = &pmo.desc_cond_lead[ (idx_C + (iprobe-1) * ct.num_blocks) * DLEN]; 
                descb = &pmo.desc_lead_cond[ (idx_C + (iprobe-1) * ct.num_blocks) * DLEN]; 
                numst = lcr[iprobe].num_states;
                numstC = ct.block_dim[idx_C];


                PDTRAN(&numst, &numstC, &one, lcr[iprobe].SCL, &ione, &ione, desca, 
                        &zero, S10, &ione, &ione, descb); 
                PDTRAN(&numst, &numstC, &one, lcr[iprobe].HCL, &ione, &ione, desca, 
                        &zero, H10, &ione, &ione, descb); 
                idx = pmo.mxllda_lead[iprobe -1] * pmo.mxlocc_cond[idx_C];
                for (i = 0; i < idx; i++)
                {
                    ch10[i] = ene * S10[i] - Ha_eV * H10[i];
                }

                Sigma_p (sigma, ch0, ch01, ch10, g, iprobe);

                for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                {
                    sigma_all[sigma_idx[iprobe - 1] + i] = sigma[i];
                }

                if(iprobe == iprobe1)
                {
                    for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                    {
                        Gamma1[i] = -2.0 * cimag(sigma[i]);
                    }
                }

                if(iprobe == iprobe2)
                {
                    for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                    {
                        Gamma2[i] = -2.0 * cimag(sigma[i]);
                    }
                }

            }  /*  end for iprobe */

            /* Construct H = ES - H */
            for (i = 0; i < ntot; i++)
            {
                H_tri[i] = creal(ene) * lcr[0].Stri[i] - lcr[0].Htri[i] * Ha_eV;
            }



            Sgreen_cond_p (H_tri, sigma_all, sigma_idx, green_C, nC, iprobe1, iprobe2);

/*   green_C store the Green function block (n2,n1)  */

            desca = &pmo.desc_cond[( n1 + n1 * ct.num_blocks) * DLEN];   /* nC_1 * nC_1 matrix */
            descb = &pmo.desc_cond[( n2 + n2 * ct.num_blocks) * DLEN];   /* nC_2 * nC_2 matrix */
            descc = &pmo.desc_cond[( n2 + n1 * ct.num_blocks) * DLEN];   /* nC_2 * nC_1 matrix */
            descd = &pmo.desc_cond[( n1 + n2 * ct.num_blocks) * DLEN];   /* nC_1 * nC_2 matrix */



            /* Gamma(C_1, C_1) * G(C_1, C_2) = temp(C_1, C_2) */
            PZGEMM (&fcd_n, &fcd_n, &nC_2, &nC_1, &nC_1, &alpha, green_C, &ione, &ione, descc,
                    Gamma1, &ione, &ione, desca, &beta, temp_matrix1, &ione, &ione, descc);

            /* temp(C_1, C_2) * Gamma(C_2, C_2) = temp2(C_1, C_2) */
            PZGEMM (&fcd_n, &fcd_c, &nC_2, &nC_2, &nC_1, &alpha, temp_matrix1, &ione, &ione, descc,
                    green_C, &ione, &ione, descc, &beta, temp_matrix2, &ione, &ione, descb);

            /* temp2(C_1, C_2) * G(C_2, C_1) = temp(C_1, C_1) */
            PZGEMM (&fcd_n, &fcd_n, &nC_2, &nC_2, &nC_2, &alpha, temp_matrix2, &ione, &ione, descb,
                    Gamma2, &ione, &ione, descb, &beta, temp_matrix1, &ione, &ione, descb);

            /* desca= &pmo.desc_lead[0]; */
            cond[iene] = pmo_trace(temp_matrix1, descb);

            ener1[iene] = creal(ene);


            /* printf (" condcond eneR, G= %f %f \n ", eneR, cond[iene]); */


        }                           /*  end for iene */

        my_barrier ();
        iene = E_POINTS;
        global_sums (cond, &iene, pct.grid_comm);

        if (pct.gridpe == 0)
        {
            sprintf(newname, "%s%s%d%d%s", ct.basename,".cond_", iprobe1, iprobe2, ".dat");
            file = fopen (newname, "w");
            for (iene = 0; iene < E_POINTS; iene++)
                fprintf (file, " %f %22.12e\n", ener1[iene], cond[iene]);
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


