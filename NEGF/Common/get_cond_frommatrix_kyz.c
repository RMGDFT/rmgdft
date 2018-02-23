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

void kpoints(int *nkp, double *kvecx, double *kvecy, double *kvecz, int *nkp_tot, double *kweight);

void get_cond_frommatrix_kyz ()
{
	int iprobe, iprobe1, iprobe2, iter;
	int iene, icond, count, insert_num;
	complex double *H_tri,*G_tri, *g;
	complex double *green_C;
	complex double *temp_matrix1, *temp_matrix2;
	complex double *Gamma1, *Gamma2, *sigma;
	double E_imag, KT, current;
	double *ener1, *cond, *ener1_temp, *cond_temp;
        int *have_cond, *have_cond_temp;
	int ntot, ndim, nC, idx_C, *sigma_idx;
	double cons, EF1, EF2, f1, f2;

	complex double *work;
	complex double ene, ctem;

	complex double alpha, beta;
	complex double one, zero;
	int i, j, idx, E_POINTS, nkp[3];
	char fcd_n = 'N', fcd_c = 'C', newname[100];
	FILE *file;
	int ione =1, *desca, *descb, *descc, *descd;
	int n1, n2, nC_1, nC_2, nC_11, nC_22, nC_max;
	int idx1, idx2;
	int num_offdiag_yz;
	int numst, numstC;
	int kp, nkp_tot;
    double cond_value, de, emin, emax;
	double *kvecx, *kvecy, *kvecz, *kweight;
    int *energy_insert_index;
    int tot_energy_point;
    int simpson_depth, simpson_loop;
    double simpson_tol, max_tol;

	/*=============== Reading input and then print them ==============*/ 

	read_cond_input (&emin, &emax, &E_POINTS, &E_imag, &KT, nkp);

    simpson_depth = ct.simpson_depth;
    simpson_tol = ct.simpson_tol;

    // if energy points are adapted by simpson's integration rule, the
    // initial E_POINTS must be 4*n+1 with n >=1 at the current code 
    //
    if(cei.energy_point_insert == 1 && simpson_depth > 0) 
    { 
        E_POINTS = (E_POINTS +2)/4 *4 +1;
    }
    
	de = (emax - emin) / (E_POINTS - 1);
    int EP;

	//  set number of kpoint in x direction be 1, 
	nkp[0] = 1;

//    if(cei.num_probe > 2 ) nkp[1] = 1;
	ntot = nkp[0] * nkp[1] * nkp[2];
	printf("\n nkp  %d %d %d", nkp[0], nkp[1], nkp[2]);

	if (ntot == 0 ) error_handler ("wrong number of kpoints in cond.in");

	my_malloc( kvecx, ntot, double );
	my_malloc( kvecy, ntot, double );
	my_malloc( kvecz, ntot, double );
	my_malloc( kweight, ntot, double );

	kpoints(nkp, kvecx, kvecy, kvecz, &nkp_tot, kweight);

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
		printf ("kpoint in x,y,z = %d %d %d\n", nkp[0], nkp[1], nkp[2]);
		printf (" total number of kpoints = %d", nkp_tot);
		printf (" kx       ky      kz      kweight" );
		for ( i = 0; i < nkp_tot; i++)
		{
			printf("\n %f, %f,  %f, %f  ", kvecx[i], kvecy[i], kvecz[i], kweight[i]);
		}
	}


   // if(nkp_tot == 1) 
    //{
     //   get_cond_frommatrix();
      //  return;
   // }

	/*======================== Reading Matrices ===============================*/


	my_malloc_init( H_tri, pmo.ntot_low, complex double );
	my_malloc_init( G_tri, pmo.ntot_low, complex double );
	my_malloc_init( lcr[0].Htri, pmo.ntot, double );
	my_malloc_init( lcr[0].Stri, pmo.ntot, double );


    nC_max = 0;
    for(i = 0; i < ct.num_blocks; i++)
    {
        idx = pmo.mxllda_cond[i] * pmo.mxlocc_cond[i]; 
        nC_max = rmg_max(nC_max, idx);
    }
    
    my_malloc_init( Gamma1, nC_max, complex double ); 
    my_malloc_init( Gamma2, nC_max, complex double ); 
    my_malloc_init( green_C, nC_max, complex double); 
    my_malloc_init( temp_matrix1, nC_max, complex double );
    my_malloc_init( temp_matrix2, nC_max, complex double );


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


    my_malloc_init(work,    12*idx, complex double );

    /*===================================================================*/

    my_malloc_init( ener1, E_POINTS * (1<<simpson_depth), double );
    my_malloc_init( cond, E_POINTS * (1<<simpson_depth), double );
    my_malloc_init( ener1_temp,  E_POINTS * (1<<simpson_depth), double );
    my_malloc_init( cond_temp, E_POINTS * (1<<simpson_depth), double );
    my_malloc_init( energy_insert_index, E_POINTS * (1<<simpson_depth), int );

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

        iprobe2 = ct.cond_probe2[icond];
        n2 = cei.probe_in_block[iprobe2 - 1];
        nC_2 = ct.block_dim[n2];


        /*===================================================================*/
        // insert another for loop to iteratively insert fine energy points between steep neighbour points. 
        // hope to get converged integration over transmission curve
        for (iene = 0; iene < E_POINTS; iene++)
        {
            ener1[iene] = emin + iene * de;
            ener1_temp[iene] = emin + iene * de;
        }

        for (iene = 0; iene < E_POINTS * (1<<simpson_depth); iene++)
            energy_insert_index[iene] = iene;

        EP = E_POINTS;
        tot_energy_point = 0;
        simpson_loop = 0;
        while(EP)  // continue while loop until EP == 0, which means there is no more points to be calculated!
        {
            for (iene =0; iene < EP; iene++)
                cond_temp[iene] = 0.0;
            for(kp = 0; kp < nkp_tot; kp++)
            {
                for (iene = pmo.myblacs; iene < EP; iene += pmo.npe_energy)
                {
                    ene = ener1_temp[iene] + I * E_imag;


                    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
                    {

                        sigma_one_energy_point(sigma, iprobe, ene, kvecy[kp], kvecz[kp], work); 


                        for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                        {
                            sigma_all[sigma_idx[iprobe - 1] + i] = sigma[i];
                        }

                        if(iprobe == iprobe1)
                        {


                            desca = &pmo.desc_cond[( n1 + n1 * ct.num_blocks) * DLEN];   /* nC_1 * nC_1 matrix */
                            pztranc(&nC_1, &nC_1, &one, sigma, &ione, &ione, desca,
                                    &zero, Gamma1, &ione, &ione, desca);
                            for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                            {
                                Gamma1[i] = I * (sigma[i] - Gamma1[i]);
                            }
                        }

                        if(iprobe == iprobe2)
                        {
                            desca = &pmo.desc_cond[( n2 + n2 * ct.num_blocks) * DLEN];   /* nC_2 * nC_2 matrix */
                            pztranc(&nC_2, &nC_2, &one, sigma, &ione, &ione, desca,
                                    &zero, Gamma2, &ione, &ione, desca);
                            for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                            {
                                Gamma2[i] = I * (sigma[i] - Gamma2[i]);
                            }
                        }

                    }  /*  end for iprobe */

                    /* Construct H = ES - H */


                    ene = ener1_temp[iene];
                    ene = ener1_temp[iene] + I * E_imag;
                    matrix_kpoint_center(H_tri, lcr[0].Stri, lcr[0].Htri, ene, kvecy[kp], kvecz[kp]);


                    Sgreen_cond_p (H_tri, G_tri, sigma_all, sigma_idx, green_C, nC, iprobe1, iprobe2);



                    desca = &pmo.desc_cond[( n1 + n1 * ct.num_blocks) * DLEN];   /* nC_1 * nC_1 matrix */
                    descb = &pmo.desc_cond[( n2 + n2 * ct.num_blocks) * DLEN];   /* nC_2 * nC_2 matrix */
                    descc = &pmo.desc_cond[( n2 + n1 * ct.num_blocks) * DLEN];   /* nC_2 * nC_1 matrix */
                    descd = &pmo.desc_cond[( n1 + n2 * ct.num_blocks) * DLEN];   /* nC_1 * nC_2 matrix */



                    /* Gamma(C_1, C_1) * G(C_1, C_2) = temp(C_1, C_2) */
                    pzgemm (&fcd_n, &fcd_n, &nC_2, &nC_1, &nC_1, &alpha, green_C, &ione, &ione, descc,
                            Gamma1, &ione, &ione, desca, &beta, temp_matrix1, &ione, &ione, descc);

                    /* temp(C_1, C_2) * Gamma(C_2, C_2) = temp2(C_1, C_2) */
                    pzgemm (&fcd_n, &fcd_c, &nC_2, &nC_2, &nC_1, &alpha, temp_matrix1, &ione, &ione, descc,
                            green_C, &ione, &ione, descc, &beta, temp_matrix2, &ione, &ione, descb);

                    /* temp2(C_1, C_2) * G(C_2, C_1) = temp(C_1, C_1) */
                    pzgemm (&fcd_n, &fcd_n, &nC_2, &nC_2, &nC_2, &alpha, temp_matrix2, &ione, &ione, descb,
                            Gamma2, &ione, &ione, descb, &beta, temp_matrix1, &ione, &ione, descb);

                    /* desca= &pmo.desc_lead[0]; */
                    cond_temp[iene] += pmo_trace(temp_matrix1, descb) * kweight[kp];

                    /* printf (" condcond eneR, G= %f %f \n ", eneR, cond[iene]); */

                }                     
            } /*  end for iene */


            my_barrier ();

            global_sums (cond_temp, &EP, pct.grid_comm);

            tot_energy_point += EP;
            for(i=0; i< EP; i++)
            {
                j = energy_insert_index[i];
                ener1[j] = ener1_temp[i];
                cond[j] = cond_temp[i];

            }
            simpson_loop++;
            if(simpson_loop > simpson_depth) 
            {
                printf("\n with Simpson depth of %d, tol = %e ", simpson_depth, max_tol);
                break;
            }

            switch(cei.energy_point_insert)
            {
                case 0:      
                    EP = 0;
                    break;
                case 1:  // Simpson integration
                    max_tol = find_new_energy_point(cond, ener1, tot_energy_point, simpson_tol, &EP, 
                            energy_insert_index, ener1_temp);
                    break;
                case 2:  // catch sharp peaks
                    find_new_energy_point_sharp_peak(cond, ener1, tot_energy_point, simpson_tol, &EP,
                            energy_insert_index, ener1_temp);
                    break;
            }
            dprintf("\n cei.energ  %d", EP);

            if(EP == 0) break;

        } //end of while loop

        if (pct.gridpe == 0)
        {
            sprintf(newname, "%s%s%d%d%s", ct.basename,".cond_", iprobe1, iprobe2, ".dat");
            file = fopen (newname, "w");

            for (iene = 0; iene < tot_energy_point; iene++)
            {
                fprintf (file, " %f %22.12e\n", ener1[iene], cond[iene]);
            }

            fclose (file);
        }
        /* find peaks of conductance, max 100 peaks, */
        int peaki=0; 
        peakNum=0; 
        for (iene = 0; iene < tot_energy_point; iene++)
        {
            if  (iene<tot_energy_point-2)
	    if  (cond[iene]<cond[iene+1] && cond[iene+2]<cond[iene+1] )
            {
                if (peaki<100)/* max 100 peaks, */
                {
                    peaks[peaki]=ener1[iene+1] ;
                    peaki++;
                }
                peakNum++;
             }
        }
        if (peakNum==0)
        {
           peaks[0]=0;
           peaks[1]=-0.5;
           peaks[2]=0.5;
           peakNum=3;
           if (pct.gridpe == 0)   printf ("\n no peak in conductance, add 3 energy points: 0, 0.5, -0.5\n");
        }
        if (pct.gridpe == 0)
        {
            printf ("\n number of peaks: %d\n", peakNum);
            for(peaki=0;peaki<peakNum;peaki++)
            printf ("\n peaks[%d]: %f\n",peaki,peaks[peaki]);
        }


        /* calculating the current */

        current = 0.0;
        EF1 = lcr[iprobe1].bias;
        EF2 = lcr[iprobe2].bias;

	if (pct.gridpe == 0)
	{
	    printf("\n ********");

	    file = fopen ("local_current.dat", "w");
	    for (iene = 1; iene < E_POINTS; iene++)
	    {
		f1 = 1.0 / (1.0 + exp ((EF1 - ener1[iene]) / KT));
		f2 = 1.0 / (1.0 + exp ((EF2 - ener1[iene]) / KT));
		current += (f2 - f1) * cond[iene] * (ener1[iene] - ener1[iene - 1]);
		fprintf(file, "\n  %f  %e", ener1[iene], (f2-f1) * cond[iene]);
	    }

	    fprintf(file, "\n  &");
	    fclose(file);

	}
	current *= 2.0 * e_C * e_C / h_SI;

        if (pct.gridpe == 0)
        {
            printf ("\n bias = %f eV    current = %e microA\n", EF1 - EF2, current * 1e6);
        }




    } /* ends of icond (num_cond_curve) loop */

    my_free(ener1);
    my_free(cond);
    my_free(ener1_temp);
    my_free(cond_temp);

    my_free(work);
    my_free(sigma);
    my_free(sigma_all);
    my_free(sigma_idx);

    my_free(Gamma1);
    my_free(Gamma2);
    my_free(temp_matrix1);
    my_free(temp_matrix2);
    my_free(green_C);
    my_free(H_tri);
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


