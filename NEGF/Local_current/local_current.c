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
#include "pmo.h"

double pmo_trace(complex double*, int*);

void kpoints(int *nkp, double *kvecx, double *kvecy, double *kvecz, int *nkp_tot, double *kweight);

void local_current ()
{
    int iprobe, iprobe1, iprobe2;
    int iene, icond;
    complex double *H_tri,*S_tri, *tot, *tott, *g;
    complex double *green_C;
    complex double *temp_matrix1, *temp_matrix2;
    complex double *Gamma1, *Gamma2, *sigma;
    REAL de, emin, emax, E_imag, KT, current;
    REAL *ener1, *cond;
    int ntot, ndim, nC, idx_C, *sigma_idx;
    REAL cons, EF1, EF2, f1, f2;

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
    double *kvecx, *kvecy, *kvecz, *kweight;
    complex double *rho_mn, *green_C_non, *Gless, *Glesst;
    int maxrow, maxcol, totrow;

/*=============== Reading input and then print them ==============*/ 

    read_cond_input (&emin, &emax, &E_POINTS, &E_imag, &KT, nkp);
    de = (emax - emin) / (E_POINTS - 1);


//  set number of kpoint in x direction be 1, 
    nkp[0] = 1;

    ntot = nkp[0] * nkp[1] * nkp[2];
    printf("\n nkp  %d %d %d", nkp[0], nkp[1], nkp[2]);

    if (ntot == 0 ) error_handler ("wrong number of kpoints in cond.in");

    my_malloc( kvecx, ntot, double );
    my_malloc( kvecy, ntot, double );
    my_malloc( kvecz, ntot, double );
    my_malloc( kweight, ntot, double );

    kpoints(nkp, kvecx, kvecy, kvecz, &nkp_tot, kweight);

    if(ct.num_blocks != 3) error_handler("local current only works for 3blocks %d", ct.num_blocks);

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
    my_malloc_init( lcr[0].Htri, ntot, REAL );
    my_malloc_init( lcr[0].Stri, ntot, REAL );

    my_malloc_init( lcr[0].Htri_yz, num_offdiag_yz *ntot, REAL );
    my_malloc_init( lcr[0].Stri_yz, num_offdiag_yz *ntot, REAL );

    maxrow = 0;
    maxcol = 0;
    totrow = 0;
    for(i = 0; i < ct.num_blocks; i++)
    {
        maxrow = max(maxrow, pmo.mxllda_cond[i]);
        maxcol = max(maxcol, pmo.mxlocc_cond[i]);
        totrow += pmo.mxllda_cond[i];
    }

    my_malloc_init( green_C_non, maxcol * totrow, complex double ); 
    my_malloc_init( rho_mn, ntot, complex double );

    i = pmo.mxllda_cond[1] * pmo.mxlocc_cond[1]; 
    my_malloc_init( Gless, i, complex double );
    my_malloc_init( Glesst, i, complex double );

    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {	
        idx = pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];
        my_malloc_init( lcr[iprobe].H00, idx, REAL );
        my_malloc_init( lcr[iprobe].S00, idx, REAL );
        my_malloc_init( lcr[iprobe].H01, idx, REAL );
        my_malloc_init( lcr[iprobe].S01, idx, REAL );

        my_malloc_init( lcr[iprobe].H00_yz, num_offdiag_yz * idx, REAL );
        my_malloc_init( lcr[iprobe].S00_yz, num_offdiag_yz * idx, REAL );
        my_malloc_init( lcr[iprobe].H01_yz, num_offdiag_yz * idx, REAL );
        my_malloc_init( lcr[iprobe].S01_yz, num_offdiag_yz * idx, REAL );
    }

    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {	
        i = cei.probe_in_block[iprobe - 1];
        idx = pmo.mxllda_cond[i] * pmo.mxlocc_lead[iprobe-1];
        my_malloc_init( lcr[iprobe].HCL, idx, REAL );
        my_malloc_init( lcr[iprobe].SCL, idx, REAL );

        my_malloc_init( lcr[iprobe].HCL_yz, num_offdiag_yz *idx, REAL );
        my_malloc_init( lcr[iprobe].SCL_yz, num_offdiag_yz *idx, REAL );
    }

    read_matrix_pp();

    if(cei.num_probe >2) 
        error_handler("cannot have ky and kz for more than two probes now");

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
    /* for left lead, the orbital index is just same as the first block in  input*/
    iprobe = 1;
    split_matrix_lead (iprobe);
    /* for right lead, the orbital index is just same as the last block in  input*/
    iprobe = 2;
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

    my_malloc_init( ener1, E_POINTS, REAL );
    my_malloc_init( cond, E_POINTS, REAL );

    alpha = 1.0;
    beta = 0.0;
    one = 1.0; 
    zero = 0.0;

    nC = ct.num_states;

    /*===================================================================*/
    file = fopen ("local_current.dat", "a");

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


        for (iene = 0; iene < E_POINTS; iene++)
        {
            ener1[iene] = emin + iene * de;
            cond[iene] = 0.0;
        }


        for (iene = pmo.myblacs; iene < E_POINTS; iene += pmo.npe_energy)
        {

            EF1 = lcr[iprobe1].bias;
            EF2 = lcr[iprobe2].bias;
            f1 = 1.0 / (1.0 + exp ((EF1 - ener1[iene]) / KT));
            f2 = 1.0 / (1.0 + exp ((EF2 - ener1[iene]) / KT));

            ene = emin + iene * de + I * E_imag;

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


                /*  calculating one row/columm of the Green function to
                 *  be used in G * Gamma * G rho_munu_p.c*/

                Sgreen_onerow (H_tri, sigma_all, sigma_idx, green_C_non, nC, iprobe1);
                rho_munu_p (rho_mn, green_C_non, Gamma1, iprobe1); 

                for (i = 0; i < pmo.mxllda_cond[1] * pmo.mxlocc_cond[1]; i++)
                {
                    Gless[i] = I * (f1-f2) * rho_mn[pmo.diag_begin[1]+i];
                }

                Sgreen_onerow (H_tri, sigma_all, sigma_idx, green_C_non, nC, iprobe2);
                rho_munu_p (rho_mn, green_C_non, Gamma2, iprobe2); 

                for (i = 0; i < pmo.mxllda_cond[1] * pmo.mxlocc_cond[1]; i++)
                {
//                    Gless[i] += I * f2 * rho_mn[pmo.diag_begin[1]+i];
                }

                desca = &pmo.desc_cond[( 1 + 1 * ct.num_blocks) * DLEN];   /* nC_1 * nC_1 matrix */
                idx1 = ct.block_dim[1];

                PZTRANU(&idx1, &idx1, &one, Gless, &ione, &ione, desca,
                        &zero, Glesst, &ione, &ione, desca);
                int *ipvc;

                my_malloc( ipvc, idx1, int );
                complex double *H, *S;
                /*
                   matrix_kpoint(ntot, H_tri, lcr[0].Htri_yz, kvecy[kp], kvecz[kp]);
                   S = &S_tri[pmo.diag_begin[1]];
                   H = &H_tri[pmo.diag_begin[1]];
                   get_inverse_block_p(S, rho_mn, ipvc, desca);

                   PZGEMM (&fcd_n, &fcd_n, &idx1, &idx1, &idx1, &alpha, H, &ione, &ione, desca,
                   rho_mn, &ione, &ione, desca, &beta, temp_matrix1, &ione, &ione, desca);
                   PZGEMM (&fcd_n, &fcd_n, &idx1, &idx1, &idx1, &alpha, rho_mn, &ione, &ione, desca,
                   H, &ione, &ione, desca, &beta, temp_matrix2, &ione, &ione, desca);

*/
                double tem;
                double tem1;
                tem = 0.0;
                for (idx1 = 128; idx1 < 134; idx1++)
                {
                    for (idx2 = 0; idx2 < pmo.mxlocc_cond[1]; idx2++)
                    {
                        i = idx2 * pmo.mxllda_cond[1] + idx1;
                        j = idx1 * pmo.mxllda_cond[1] + idx2;


                        tem += creal( (Glesst[i]- Gless[i]) * H_tri[pmo.diag_begin[1] +i] );
                        tem1 = creal( (Glesst[i]- Gless[i]) * H_tri[pmo.diag_begin[1] +i] );
                        printf("\n %f  %e  %d ", creal(ene), tem1, idx2);
                        // tem += creal( (Gless[i]- Glesst[i]) * H[i] * Ha_eV);
                        // tem += creal( Gless[i] *temp_matrix2[i] - Gless[j] * temp_matrix1[j]) * Ha_eV;
                    }

                }

                fprintf(file, "\n  %f  %e", creal(ene), tem);

            }                           /*  end for iene */
        }

        my_barrier ();


    } /* ends of icond (num_cond_curve) loop */
    fclose(file);

    my_barrier ();

}

