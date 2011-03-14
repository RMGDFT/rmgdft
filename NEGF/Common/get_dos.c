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


void get_dos (STATE * states)
{
    int iprobe;
    int iene, st1;
    REAL eneR, eneI;
    doublecomplex *tot, *tott, *g;
    REAL *green_tem, *green_C;
    doublecomplex *sigma, *sigma_all;
    REAL *temp_matrix_tri, *temp_matrix, *matrix_product;
    REAL de, emin, emax;

    int nC;
    int i, j, *sigma_idx, idx_C;
    char fcd_n = 'N', fcd_c = 'C';
    FILE *file;

    int ntot, ndim;
    int ii, jj, kk, xoff, yoff, zoff;
    REAL *Green_store, *rho_energy, *rho_energy2;
    int root_pe, idx, ix, iy, iz;

    int E_POINTS, kpoint;
    double E_imag, KT;
    int FPYZ = FPY0_GRID * FPZ0_GRID;
    int nx1, nx2, ny1, ny2, nz1, nz2; 

    read_cond_input (&emin, &emax, &E_POINTS, &E_imag, &KT, &kpoint);
    de = (emax - emin) / (E_POINTS - 1);



    /* store the imaginary part of Green function Gc on each processor
     * then calculate the energy dependent charge density 
     * and stored in rho_energy after average in the yz plane
     */

    nC = ct.num_states;

/*========================== Reading Matrices =======================*/

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


/*======================== Allocate memory for sigma =================*/

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


    my_malloc_init( green_tem, 2 * idx, REAL );
    my_malloc_init( green_C, 2 * ntot, REAL );
    /*st1 = (E_POINTS + NPES - 1) / NPES;*/
    st1 = (E_POINTS + pmo.npe_energy - 1) / pmo.npe_energy;
   
    my_malloc_init( Green_store, st1 * ntot, REAL );

/*===================================================================*/

    nx1 = cei.dos_window_start[0] * RHO_NX;
    nx2 = cei.dos_window_end[0] * RHO_NX;
    ny1 = cei.dos_window_start[1] * RHO_NY;
    ny2 = cei.dos_window_end[1] * RHO_NY;
    nz1 = cei.dos_window_start[2] * RHO_NZ;
    nz2 = cei.dos_window_end[2] * RHO_NZ;
                                                                                              
                                                                                              
                                                                                              
                                                                                              
    my_malloc_init( rho_energy, E_POINTS * FNX_GRID, REAL );
    if (cei.num_probe > 2)
        my_malloc_init( rho_energy2, E_POINTS * FNY_GRID, REAL );
                                                                                              
                                                                                              
                                                                                              
    pe2xyz (pct.gridpe, &ii, &jj, &kk);
    xoff = ii * FPX0_GRID;
    yoff = jj * FPY0_GRID;
    zoff = kk * FPZ0_GRID;


/*===================================================================*/



    idx = 0;
    /*for (iene = pct.gridpe; iene < E_POINTS; iene += NPES)*/
    for (iene = pmo.myblacs; iene < E_POINTS; iene += pmo.npe_energy)
    {
        eneR = emin + iene * de;
        eneI = 0.0005;
        printf ("\n energy point %d %f\n", iene, eneR);


        /* sigma is a complex matrix with dimension ct.num_states * ct.num_states 
         * it sums over all probes
         * tot, tott, green_tem is also a complex matrix, 
         * their memory should be the maximum of probe dimensions, lcr[1,...].num_states
         */
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


        }                       /*  end for iprobe */




        /*Sgreen_c_wang (lcr[0].Htri, lcr[0].Stri, sigma_all, sigma_idx, 
                      eneR, eneI, (doublecomplex *) green_C, nC);*/
        Sgreen_c_p (lcr[0].Htri, lcr[0].Stri, sigma_all, sigma_idx, 
                      eneR, eneI, (doublecomplex *) green_C); 


        for (st1 = 0; st1 < ntot; st1++)
        {
            Green_store[idx + st1] = -green_C[st1 * 2 + 1];
            /*printf (" checkit  = %d %d %f \n", iene, idx + st1, Green_store[idx + st1]);*/
        }
        idx += ntot;

    }                           /*  end for iene */




/*===================================================================*/


    /*for (iene = pmo.myblacs; iene < E_POINTS; iene += pmo.npe_energy)*/
    for (iene = 0; iene < E_POINTS; iene++)
    {
            printf ("hello .... %d\n", iene);
                                                                                                                  
        root_pe = iene % pmo.npe_energy;
        idx = iene / pmo.npe_energy;


        for (st1 = 0; st1 < ntot; st1++)
        {
            lcr[0].density_matrix_tri[st1] = Green_store[idx * ntot + st1];
        }
    
                                                                                                             
        idx = ntot;
        MPI_Bcast (lcr[0].density_matrix_tri, idx, MPI_DOUBLE, root_pe,
COMM_EN1);
                                                                                                                
        get_new_rho_soft (states, rho);
                                                                                                                  

        for (ix = 0; ix < FPX0_GRID; ix++)
        {
            for (iy = 0; iy < FPY0_GRID; iy++)
            {
                if((iy + yoff) >= ny1 && (iy + yoff) <= ny2)
                {
                    for (iz = 0; iz < FPZ0_GRID; iz++)
                    {
                        if((iz + zoff) >= nz1 && (iz + zoff) <= nz2)
                        {
                            idx = iz + iy * FPZ0_GRID + ix * FPYZ;
                            rho_energy[iene * FNX_GRID + ix + xoff] += rho[idx];
                        }
                    }
                }
            }
        }
                                                                                                                  
                                                                                                                  
        if (cei.num_probe > 2)
        {
            for (iy = 0; iy < FPY0_GRID; iy++)
            {
                for (ix = 0; ix < FPX0_GRID; ix++)
                {
                    if((ix + xoff) >= nx1 && (ix + xoff) <= nx2)
                    {
                        for (iz = 0; iz < FPZ0_GRID; iz++)
                        {
                            if((iz + zoff) >= nz1 && (iz + zoff) <= nz2)
                            {
                                idx = iz + iy * FPZ0_GRID + ix * FPYZ;
                                rho_energy2[iene * FNY_GRID + iy + yoff] += rho[idx];
                            }
                        }
                    }
                }
            }
        }
                                                                                                                  
        my_barrier ();
    }



    iene = E_POINTS * FNX_GRID;
    global_sums (rho_energy, &iene);
    if (pct.gridpe == 0)
    {
        double dx = ct.celldm[0] / NX_GRID;
        double x0 = 0.5 * ct.celldm[0];
                                                                                                                  
        file = fopen ("dos.dat", "w");
        fprintf (file, "#     x[a0]      E[eV]          dos\n\n");
        for (iene = 0; iene < E_POINTS; iene++)
        {
            eneR = emin + iene * de;
                                                                                                                  
            for (ix = 0; ix < NX_GRID; ix++)
            {
                                                                                                                  
                fprintf (file, " %10.6f %10.6f %12.6e\n",
                         ix * dx - x0, eneR, rho_energy[iene * FNX_GRID + ix * RHO_NX]);
            }
            fprintf (file, "\n");
        }
                                                                                                                  
        fclose (file);
    }

    if (cei.num_probe > 2)
    {
                                                                                                                  
        iene = E_POINTS * FNY_GRID;
        global_sums (rho_energy2, &iene);
        if (pct.gridpe == 0)
        {
            double y = ct.celldm[1] * ct.celldm[0];
            double dy = y / NY_GRID;
            double y0 = 0.5 * y;
                                                                                                                  
            file = fopen ("dos2.dat", "w");
            fprintf (file, "#     y[b0]      E[eV]          dos\n\n");
            for (iene = 0; iene < E_POINTS; iene++)
            {
                eneR = emin + iene * de;
                                                                                                                  
                for (iy = 0; iy < NY_GRID; iy++)
                {
                                                                                                                  
                    fprintf (file, " %10.6f %10.6f %12.6e\n",
                            iy * dy - y0, eneR, rho_energy2[iene * FNY_GRID + iy * RHO_NY]);
                }
                fprintf (file, "\n");
            }
                                                                                                                  
            fclose (file);
        }
        my_free(rho_energy2);
    }
                                                                                                                  
                                                                                                                  
    my_barrier ();
    fflush (NULL);

/*===============================*/
    my_free(tot);
    my_free(tott);
    my_free(g);
    my_free(green_tem);
    my_free(sigma_all);
    my_free(sigma_idx);
    my_free(green_C);

    my_free(lcr[0].Htri);
    my_free(lcr[0].Stri);
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        my_free(lcr[iprobe].H00);
        my_free(lcr[iprobe].S00);
        my_free(lcr[iprobe].H01);
        my_free(lcr[iprobe].S01);
    }

/*===============================*/

    my_free(rho_energy); 
    my_free(Green_store);



}
