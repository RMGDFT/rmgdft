/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "main.h"
#include "init_var_negf.h"
#include "LCR.h"

#define 	LDEBUG 	0

void multi_GHG_munu (rmg_double_t *GHG_tri, rmg_double_t *GHG_en_tri)
{
    int iprobe, jprobe, iene;
    int st1, st2, idx_sigma;
    int nC, nL, i, ntot, ion;
    int idx1, idx2, nmax, size;
    complex double ene, weight;
    rmg_double_t wmn, wmn1, wmn2;
    complex double *green_C, *complex_H, *temp_matrix_tri, *temp_matrix1;
    complex double *sigma_sum1, *sigma_sum2, *sigma_L;
    complex double one, zero;
    double time1, time2, time3, time4;

    time1 = my_crtc ();

    one = 1.0;
    zero = 0.0;

    /*  allocate memory for sigma_sum1, only the first block. others's value is zero, don't need to save it*/
    nL = lcr[1].num_states;
    my_malloc_init( sigma_sum1, nL * nL, complex double );
    if (nL != ct.block_dim[0])
    {
        printf (" lcr[1].num_states & ct.block_dim[0] are unequal \n");
    }

    /*  allocate memory for sigma_sum2, only the last block. others's value is zero, don't need to save it*/
    nL = lcr[2].num_states;
    my_malloc_init( sigma_sum2, nL * nL, complex double );
    if (nL != ct.block_dim[ct.num_blocks - 1])
    {
        printf (" lcr[2].num_states & ct.block_dim[%d] are unequal \n", ct.num_blocks - 1);
    }

    /*  allocate memory for grenn_C_tri is tri-diagonal */
    ntot = 0;
    for (i = 0; i < ct.num_blocks - 1; i++)
        ntot += (ct.block_dim[i] * ct.block_dim[i] + ct.block_dim[i] * ct.block_dim[i + 1]);
    ntot += ct.block_dim[ct.num_blocks - 1] * ct.block_dim[ct.num_blocks - 1];
    nmax = 0;
    for (i = 0; i < ct.num_blocks; i++)
    {
        if (nmax < ct.block_dim[i])
            nmax = ct.block_dim[i];
    }

    my_malloc_init( temp_matrix_tri, ntot, complex double );

    size = ct.num_states * ct.num_states;
    my_malloc_init( green_C, size, complex double );
    my_malloc_init( complex_H, size, complex double );
    my_malloc_init( temp_matrix1, size, complex double );

    idx_sigma = 0;

    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {

        /*   parallel for the processor on energy grid */
        /*   make sure the energy points on each process are same, this guarantee 
         *   the parallelization of the calculation of par_Hij_tri and par_Sij_tri*/
        for (iene = pct.gridpe; iene < lcr[iprobe].nenergy; iene += NPES)
        {

            ene = lcr[iprobe].ene[iene];
            weight = lcr[iprobe].weight[iene];


            /* sigma is a complex matrix with dimension ct.num_states * ct.num_states 
             * it sums over all probes
             * tot, tott, green_tem is also a complex matrix, 
             * their memory should be the maximum of probe dimensions, lcr[1,...].num_states
             */


            for (jprobe = 1; jprobe <= cei.num_probe; jprobe++)
            {


                nL = lcr[jprobe].num_states;
                nC = ct.num_states;

                if (jprobe == 1)
                {
                    for (st1 = 0; st1 < nL; st1++)
                        for (st2 = 0; st2 < nL; st2++)
                        {
                            idx1 = st1 * nL + st2;
                            idx2 = st1 * nL + st2;
                            sigma_sum1[idx1] = sigma_all[idx_sigma + idx2];
                        }

                    idx_sigma += nL * nL;
                }
                else if (jprobe == 2)
                {
                    for (st1 = 0; st1 < nL; st1++)
                        for (st2 = 0; st2 < nL; st2++)
                        {
                            idx1 = st1 * nL + st2;
                            idx2 = st1 * nL + st2;
                            sigma_sum2[idx1] = sigma_all[idx_sigma + idx2];
                        }
                    idx_sigma +=  nL * nL;
                }
                else
                {
                    error_handler ("haha");
                }

            }


    	    time3 = my_crtc ();

            /* Calculating the Green function */
            Sgreen_c (lcr[0].Htri, lcr[0].Stri, sigma_sum1, sigma_sum2, ene, green_C, nC);

    	    time4 = my_crtc ();
    	    rmg_timings (NLFORCE_GREE_F, (time4 - time3));

            /* Expanding Htri to complex space and temporarily stored in temp_matrix_tri*/
            for (idx1 = 0; idx1 < ntot; idx1++)
            {
                if (ct.runflag == 113)
                {
                    wmn1 = lcr[1].lcr_ne[0].density_matrix_ne_tri[st1];
                    wmn2 = lcr[2].lcr_ne[0].density_matrix_ne_tri[st1];
                    wmn1 = wmn1 * wmn1;
                    wmn2 = wmn2 * wmn2;

                    wmn = wmn2 / (wmn1 + wmn2);
                    if (wmn1 + wmn2 < 0.0000001)
                    wmn = 0.5;

                    if (iprobe == 1) temp_matrix_tri[idx1] = wmn * lcr[0].Htri[idx1];
                    else temp_matrix_tri[idx1] = (1.0 - wmn) * lcr[0].Htri[idx1];
                }
                else
                {
                    temp_matrix_tri[idx1] = lcr[0].Htri[idx1];
                }
            }

            tri_to_whole_complex( temp_matrix_tri, complex_H, ct.num_blocks, ct.block_dim);

            /* Calculating Green_fun^T * Hamiltonin * Green_fun^T */
            ZGEMM ("T", "N", &nC, &nC, &nC, &one, green_C, &nC, complex_H, &nC, &zero, temp_matrix1, &nC);
            ZGEMM ("N", "T", &nC, &nC, &nC, &one, temp_matrix1, &nC, green_C, &nC, &zero, complex_H, &nC);

            /* Transfering (G^T * H * G^T) matrix to tridiganol matrix and stored temporarily in temp_matrix_tri */
            whole_to_tri_complex (temp_matrix_tri, complex_H, ct.num_blocks, ct.block_dim);

            /* Now performing the integration */
            for(idx1 = 0; idx1 < ntot; idx1++)
            {
                GHG_tri[idx1] +=  cimag(weight * temp_matrix_tri[idx1]) ;
                GHG_en_tri[idx1] += cimag ( weight * temp_matrix_tri[idx1] * ene);
            }


        }

        if (ct.runflag == 111 | ct.runflag == 112)
            break;

    }

    /* Computing the part from the non-equilibrium term  used for force calculation
     * Here GHG actually is G^T * H * Rho^T 
     */
    if (ct.runflag == 113)
    {
        my_malloc_init( sigma_L, size, complex double );

        for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
        {

            /*   parallel for the processor on energy grid */
            for (iene = pct.gridpe; iene < lcr[iprobe].lcr_ne[0].nenergy_ne; iene += NPES)
            {

                ene = lcr[iprobe].lcr_ne[0].ene_ne[iene];
                weight = lcr[iprobe].lcr_ne[0].weight_ne[iene];


                /* sigma is a complex matrix with dimension ct.num_states * ct.num_states 
                 * it sums over all probes
                 * tot, tott, green_tem is also a complex matrix, 
                 * their memory should be the maximum of probe dimensions, lcr[1,...].num_states
                 */


                for (jprobe = 1; jprobe <= cei.num_probe; jprobe++)
                {
                    nL = lcr[jprobe].num_states;
                    nC = ct.num_states;

                    if (jprobe == 1)
                    {
                        for (st1 = 0; st1 < nL; st1++)
                            for (st2 = 0; st2 < nL; st2++)
                            {
                                idx1 = st1 * nL + st2;
                                idx2 = st1 * nL + st2;
                                sigma_sum1[idx1] = sigma_all[idx_sigma + idx2];
                            }

                        idx_sigma +=  nL * nL;
                    }
                    else if (jprobe == 2)
                    {
                        for (st1 = 0; st1 < nL; st1++)
                            for (st2 = 0; st2 < nL; st2++)
                            {
                                idx1 = st1 * nL + st2;
                                idx2 = st1 * nL + st2;
                                sigma_sum2[idx1] = sigma_all[idx_sigma + idx2];
                            }
                        idx_sigma +=  nL * nL;
                    }
                    else
                    {
                        error_handler ("haha");
                    }



                }

                for (st1 = 0; st1 < size; st1++) 
                {
                    sigma_L[st1] = 0.0;
                }

                nL = lcr[iprobe].num_states;
                nC = ct.num_states;

                if (iprobe == 1)
                {
                    for (st1 = 0; st1 < nL; st1++)
                        for (st2 = 0; st2 < nL; st2++)
                        {
                            idx1 = st1 * nC + st2;
                            idx2 = st1 * nL + st2;
                            sigma_L[idx1] = sigma_sum1[idx2];
                        }

                }
                else if (iprobe == 2)
                {
                    for (st1 = 0; st1 < nL; st1++)
                        for (st2 = 0; st2 < nL; st2++)
                        {
                            idx1 = (nC - nL + st1) * nC + (nC - nL + st2);
                            idx2 = st1 * nL + st2;
                            sigma_L[idx1] = sigma_sum2[idx2];
                        }
                }



                /* Calculating the Green function */
                Sgreen_c (lcr[0].Htri, lcr[0].Stri, sigma_sum1, sigma_sum2, ene, green_C, nC);

                /* Expanding Htri to complex space and temporarily stored in temp_matrix_tri*/
                for (idx1 = 0; idx1 < ntot; idx1++)
                {
                    wmn1 = lcr[1].lcr_ne[0].density_matrix_ne_tri[st1];
                    wmn2 = lcr[2].lcr_ne[0].density_matrix_ne_tri[st1];
                    wmn1 = wmn1 * wmn1;
                    wmn2 = wmn2 * wmn2;

                    wmn = wmn2 / (wmn1 + wmn2);
                    if (wmn1 + wmn2 < 0.0000001)
                        wmn = 0.5;

                    if ( iprobe == 2) temp_matrix_tri[idx1] = wmn * lcr[0].Htri[idx1];
                    else temp_matrix_tri[idx1] = (1.0 - wmn) * lcr[0].Htri[idx1];
                }

                tri_to_whole_complex( temp_matrix_tri, complex_H, ct.num_blocks, ct.block_dim);

                /* Calculating Pi * G^T * H * Rho^T = G^T * H * (G**) * Sigma^T * G^T 
                 * since Rho = 1/Pi * G * Sigma * G^+     
                 */
                ZGEMM ("T", "N", &nC, &nC, &nC, &one, green_C, &nC, complex_H, &nC, &zero, temp_matrix1, &nC);
                ZGEMM ("N", "C", &nC, &nC, &nC, &one, temp_matrix1, &nC, green_C, &nC, &zero, complex_H, &nC);
                ZGEMM ("N", "T", &nC, &nC, &nC, &one, complex_H, &nC, sigma_L, &nC, &zero, temp_matrix1, &nC);
                ZGEMM ("N", "T", &nC, &nC, &nC, &one, temp_matrix1, &nC, green_C, &nC, &zero, complex_H, &nC);

                /* Transfering (G^T * H * Rho^T) matrix to tridiganol matrix and stored temporarily in temp_matrix_tri */
                whole_to_tri_complex (temp_matrix_tri, complex_H, ct.num_blocks, ct.block_dim);

                /* Now performing the integration */
                for(idx1 = 0; idx1 < ntot; idx1++)
                {
                    GHG_tri[idx1] -= 2.0 /PI * creal(weight * temp_matrix_tri[idx1]);
                    GHG_en_tri[idx1] -= 2.0 /PI * creal(weight * temp_matrix_tri[idx1] * ene); 
                }


            }

        }

        my_free(sigma_L);

    }

    global_sums (GHG_tri, &ntot, pct.grid_comm);
    global_sums (GHG_en_tri, &ntot, pct.grid_comm);


    for (idx1 = 0; idx1 < ntot; idx1++)
    {
        GHG_tri[idx1] /= PI;
        GHG_en_tri[idx1] /= PI;
    }

    my_free(sigma_sum1);
    my_free(sigma_sum2);
    my_free(temp_matrix_tri);
    my_free(green_C);
    my_free(complex_H);
    my_free(temp_matrix1);

    time2 = my_crtc ();
    rmg_timings (PAR_D_GHG, (time2 - time1));

}
