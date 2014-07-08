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
	complex double *H_tri,*S_tri, *tot, *tott, *g;
	complex double *green_C;
	complex double *temp_matrix1, *temp_matrix2;
	complex double *Gamma1, *Gamma2, *sigma;
	double E_imag, KT, current;
	double *ener1, *cond, *ener1_temp, *cond_temp;
        int *have_cond, *have_cond_temp;
	int ntot, ndim, nC, idx_C, *sigma_idx;
	double cons, EF1, EF2, f1, f2;

	complex double *ch0, *ch01, *ch10;
	complex double *H10, *S10, *H01, *H00, *S01, *S00;
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


	num_offdiag_yz = 9;


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
	my_malloc_init( S_tri, ntot, complex double );
	my_malloc_init( lcr[0].Htri, ntot, double );
	my_malloc_init( lcr[0].Stri, ntot, double );

	my_malloc_init( lcr[0].Htri_yz, num_offdiag_yz *ntot, double );
	my_malloc_init( lcr[0].Stri_yz, num_offdiag_yz *ntot, double );

	for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
	{	
		idx = pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];
		my_malloc_init( lcr[iprobe].H00, idx, double );
		my_malloc_init( lcr[iprobe].S00, idx, double );
		my_malloc_init( lcr[iprobe].H01, idx, double );
		my_malloc_init( lcr[iprobe].S01, idx, double );

		my_malloc_init( lcr[iprobe].H00_yz, num_offdiag_yz * idx, double );
		my_malloc_init( lcr[iprobe].S00_yz, num_offdiag_yz * idx, double );
		my_malloc_init( lcr[iprobe].H01_yz, num_offdiag_yz * idx, double );
		my_malloc_init( lcr[iprobe].S01_yz, num_offdiag_yz * idx, double );
	}

	for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
	{	
		i = cei.probe_in_block[iprobe - 1];
		idx = pmo.mxllda_cond[i] * pmo.mxlocc_lead[iprobe-1];
		my_malloc_init( lcr[iprobe].HCL, idx, double );
		my_malloc_init( lcr[iprobe].SCL, idx, double );

		my_malloc_init( lcr[iprobe].HCL_yz, num_offdiag_yz *idx, double );
		my_malloc_init( lcr[iprobe].SCL_yz, num_offdiag_yz *idx, double );
	}

	read_matrix_pp();

//	if(cei.num_probe >2) 
//		error_handler("cannot have ky and kz for more than two probes now");

	/*===========   all matrix elements are splitted according to ========
	 * (1) if two overlaping orbitals are in the same cell, the value stored
	 *      in Matrix_yz[0+ij],  ij is the index of the matrix element.
	 * (2) if two orbitals overlaps after one translate by unit vector in y
	 *     and/or z, the values are stored in Matrix_yz[K+ij]
	 *     K = size of matrix  * m 
	 *     m = 1 for +y translation
	 *     m = 2 for -y translation
	 *     m = 3 for +z translation
	 *     m = 4 for -z translation
	 *     m = 5 for +y+z translation
	 *     m = 6 for +y-z translation
	 *     m = 7 for -y+z translation
	 *     m = 8 for -y-z translation
	 */


	/* for center part, the orbital index is just same as input*/
	split_matrix_center ();  

	for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    	split_matrix_lead (iprobe);


	/*=================== Allocate memory for sigma ======================*/

	idx = 0;
	for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
	{
		idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
		idx = max(idx, pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]);
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
		idx = max(idx, pmo.mxllda_cond[idx_C] * pmo.mxlocc_lead[iprobe-1]);
	}

	my_malloc_init( tot,  idx, complex double );
	my_malloc_init( tott, idx, complex double );
	my_malloc_init( g,    idx, complex double );
	my_malloc_init( ch0,  idx, complex double );
	my_malloc_init( ch01,  idx, complex double );
	my_malloc_init( ch10,  idx, complex double );
	my_malloc_init( H00,  idx, complex double );
	my_malloc_init( H01,  idx, complex double );
	my_malloc_init( H10,  idx, complex double );
	my_malloc_init( S00,  idx, complex double );
	my_malloc_init( S01,  idx, complex double );
	my_malloc_init( S10,  idx, complex double );

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
        nC_11 = pmo.mxllda_cond[n1] * pmo.mxlocc_cond[n1];

        iprobe2 = ct.cond_probe2[icond];
        n2 = cei.probe_in_block[iprobe2 - 1];
        nC_2 = ct.block_dim[n2];
        nC_22 = pmo.mxllda_cond[n2] * pmo.mxlocc_cond[n2];

        my_malloc_init( Gamma1, nC_11, complex double ); 
        my_malloc_init( Gamma2, nC_22, complex double ); 

        nC_max = max(nC_11, nC_22);

        my_malloc_init( green_C, nC_max, complex double); 
        my_malloc_init( temp_matrix1, nC_max, complex double );
        my_malloc_init( temp_matrix2, nC_max, complex double );

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

        if (pct.gridpe == 0)
        {
            sprintf(newname, "%s%d%d%s", "cond_", iprobe1, iprobe2, ".dat1");
            file = fopen (newname, "w");
        }

        EP = E_POINTS;
        tot_energy_point = 0;
        simpson_loop = 0;
        while(EP)  // continue while loop until EP == 0, which means there is no more points to be calculated!
        {
            for (iene =0; iene < EP; iene++)
                cond_temp[iene] = 0.0;
            for (iene = pmo.myblacs; iene < EP; iene += pmo.npe_energy)
            {
                ene = ener1_temp[iene] + I * E_imag;

                for(kp = 0; kp < nkp_tot; kp++)
                {


                    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
                    {

                        idx = pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];

                        matrix_kpoint(idx, S01, lcr[iprobe].S01_yz, kvecy[kp], kvecz[kp]);
                        matrix_kpoint(idx, H01, lcr[iprobe].H01_yz, kvecy[kp], kvecz[kp]);
                        matrix_kpoint(idx, S00, lcr[iprobe].S00_yz, kvecy[kp], kvecz[kp]);
                        matrix_kpoint(idx, H00, lcr[iprobe].H00_yz, kvecy[kp], kvecz[kp]);

                        desca = &pmo.desc_lead[ (iprobe-1) * DLEN];

                        numst = lcr[iprobe].num_states;

                        PZTRANC(&numst, &numst, &one, S01, &ione, &ione, desca,
                                &zero, S10, &ione, &ione, desca);
                        PZTRANC(&numst, &numst, &one, H01, &ione, &ione, desca,
                                &zero, H10, &ione, &ione, desca);



                        for (i = 0; i < idx; i++)
                        {
                            ch0 [i] = ene * S00[i] - Ha_eV * H00[i];
                            ch01[i] = ene * S01[i] - Ha_eV * H01[i];
                            ch10[i] = ene * S10[i] - Ha_eV * H10[i];
                        }


                        Stransfer_p (tot, tott, ch0, ch01, ch10,iprobe);

                        Sgreen_p (tot, tott, ch0, ch01, g, iprobe);


                        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
                        idx = pmo.mxllda_cond[idx_C] * pmo.mxlocc_lead[iprobe-1];

                        matrix_kpoint(idx, S01, lcr[iprobe].SCL_yz, kvecy[kp], kvecz[kp]);
                        matrix_kpoint(idx, H01, lcr[iprobe].HCL_yz, kvecy[kp], kvecz[kp]);

                        for (i = 0; i < idx; i++)
                        {
                            ch01[i] = creal(ene) * S01[i] - Ha_eV * H01[i];
                        }

                        desca = &pmo.desc_cond_lead[ (idx_C + (iprobe-1) * ct.num_blocks) * DLEN];
                        descb = &pmo.desc_lead_cond[ (idx_C + (iprobe-1) * ct.num_blocks) * DLEN];
                        numst = lcr[iprobe].num_states;
                        numstC = ct.block_dim[idx_C];


                        PZTRANC(&numst, &numstC, &one, S01, &ione, &ione, desca,
                                &zero, S10, &ione, &ione, descb);
                        PZTRANC(&numst, &numstC, &one, H01, &ione, &ione, desca,
                                &zero, H10, &ione, &ione, descb);
                        idx = pmo.mxllda_lead[iprobe -1] * pmo.mxlocc_cond[idx_C];
                        for (i = 0; i < idx; i++)
                        {
                            ch10[i] = creal(ene) * S10[i] - Ha_eV * H10[i];
                        }

                        Sigma_p (sigma, ch0, ch01, ch10, g, iprobe);




                        for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                        {
                            sigma_all[sigma_idx[iprobe - 1] + i] = sigma[i];
                        }

                        if(iprobe == iprobe1)
                        {


                            desca = &pmo.desc_cond[( n1 + n1 * ct.num_blocks) * DLEN];   /* nC_1 * nC_1 matrix */
                            PZTRANC(&nC_1, &nC_1, &one, sigma, &ione, &ione, desca,
                                    &zero, Gamma1, &ione, &ione, desca);
                            for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                            {
                                Gamma1[i] = I * (sigma[i] - Gamma1[i]);
                            }
                        }

                        if(iprobe == iprobe2)
                        {
                            desca = &pmo.desc_cond[( n2 + n2 * ct.num_blocks) * DLEN];   /* nC_2 * nC_2 matrix */
                            PZTRANC(&nC_2, &nC_2, &one, sigma, &ione, &ione, desca,
                                    &zero, Gamma2, &ione, &ione, desca);
                            for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                            {
                                Gamma2[i] = I * (sigma[i] - Gamma2[i]);
                            }
                        }

                    }  /*  end for iprobe */

                    /* Construct H = ES - H */

                    matrix_kpoint(ntot, H_tri, lcr[0].Htri_yz, kvecy[kp], kvecz[kp]);
                    matrix_kpoint(ntot, S_tri, lcr[0].Stri_yz, kvecy[kp], kvecz[kp]);


                    for (i = 0; i < ntot; i++)
                    {
                        H_tri[i] = creal(ene) * S_tri[i] - H_tri[i] * Ha_eV;
                    }



                    Sgreen_cond_p (H_tri, sigma_all, sigma_idx, green_C, nC, iprobe1, iprobe2);





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
            for (iene = 0; iene < tot_energy_point; iene++)
            {
                fprintf (file, " %f %22.12e\n", ener1[iene], cond[iene]);
            }

            fclose (file);
        }


        /* calculating the current */

        current = 0.0;
        EF1 = lcr[iprobe1].bias;
        EF2 = lcr[iprobe2].bias;

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
    my_free(ener1_temp);
    my_free(cond_temp);

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


