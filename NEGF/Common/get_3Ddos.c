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
#include <complex.h>

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"


void get_3Ddos (STATE * states, double EMIN, double EMAX, int EPoints, int number)
{
    int iprobe, idx_e;
    int iene, st1;
    complex double *sigma;
    double de, emin, emax;

    complex double *work, *green_C;
    complex double ene;
    int kp, i, *sigma_idx, idx_C;
    FILE *file;

    int ntot, ndim;
    int  xoff, yoff, zoff;
    double *Green_store, *rho_energy, *rho_energy2;
    int root_pe, idx, ix, iy, iz;

    int E_POINTS, nkp[3];
    double E_imag, KT;
    int FPYZ = get_FPY0_GRID() * get_FPZ0_GRID();
    int nx1, nx2, ny1, ny2, nz1, nz2; 
    double *kvecx, *kvecy, *kvecz, *kweight;


    read_cond_input (&emin, &emax, &E_POINTS, &E_imag, &KT, nkp);
    de = (EMAX - EMIN) / (EPoints - 1);

    //  set number of kpoint in x direction be 1, 
    nkp[0] = 1;

    if(cei.num_probe > 2 ) nkp[1] = 1;
    int nkp_tot = nkp[0] * nkp[1] * nkp[2];
    printf("\n nkp  %d %d %d", nkp[0], nkp[1], nkp[2]);

    if (nkp_tot == 0 ) error_handler ("wrong number of kpoints in cond.in");

    my_malloc( kvecx, nkp_tot, double );
    my_malloc( kvecy, nkp_tot, double );
    my_malloc( kvecz, nkp_tot, double );
    my_malloc( kweight, nkp_tot, double );

    kpoints(nkp, kvecx, kvecy, kvecz, &nkp_tot, kweight);




    /* store the imaginary part of Green function Gc on each processor
     * then calculate the energy dependent charge density 
     * and stored in rho_energy after average in the yz plane
     */


/*========================== Reading Matrices =======================*/

    ndim = 0;
    for (i = 0; i < ct.num_blocks; i++)
    {
        ndim += ct.block_dim[i];
    }
    if (ndim != ct.num_states)
    {
        printf (" %d %d ndim not equal to nC in get_cond_frommatrix\n", ndim, ct.num_states);
        exit (0);
    }


    ntot = pmo.ntot;
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


/*======================== Allocate memory for sigma =================*/

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
                                                                                
    my_malloc_init( work,  12 * idx, complex double );


    my_malloc_init( green_C, pmo.ntot_low, complex double );

    st1 = ( E_POINTS + pmo.npe_energy-1)/pmo.npe_energy;
    my_malloc_init( Green_store, st1 * pmo.ntot, double );
    double *density_matrix;
    my_malloc_init( density_matrix, pmo.ntot, double );
/*===================================================================*/

    nx1 = cei.dos_window_start[0] * get_FG_RATIO();
    nx2 = cei.dos_window_end[0] * get_FG_RATIO();
    ny1 = cei.dos_window_start[1] * get_FG_RATIO();
    ny2 = cei.dos_window_end[1] * get_FG_RATIO();
    nz1 = cei.dos_window_start[2] * get_FG_RATIO();
    nz2 = cei.dos_window_end[2] * get_FG_RATIO();
                                                                                              
                                                                                              
                                                                                              
                                                                                              
    my_malloc_init( rho_energy, get_FNX_GRID() * get_FNY_GRID() * get_FNZ_GRID(), double );
                                                                                              
                                                                                              
                                                                                              
    xoff = get_FPX_OFFSET();
    yoff = get_FPY_OFFSET();
    zoff = get_FPZ_OFFSET();


/*===================================================================*/



    idx_e = 0;
    //for (iene = 0; iene < EPoints; iene ++)
    for (iene = pmo.myblacs; iene < EPoints; iene += pmo.npe_energy)
    {
        ene = EMIN + iene * de + I * E_imag;

        /* sigma is a complex matrix with dimension ct.num_states * ct.num_states 
         * it sums over all probes
         * tot, tott, green_tem is also a complex matrix, 
         * their memory should be the maximum of probe dimensions, lcr[1,...].num_states
         */
        for(kp = 0; kp < nkp_tot; kp++)
        {
            for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
            {

                sigma_one_energy_point(sigma, iprobe, ene, kvecy[kp], kvecz[kp], work);

                for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                {
                    sigma_all[sigma_idx[iprobe - 1] + i] = sigma[i];
                }

            }                       /*  end for iprobe */





            Sgreen_c_p (lcr[0].Htri, lcr[0].Stri, sigma_all, sigma_idx, 
                    ene, green_C); 


            for (st1 = 0; st1 < ntot; st1++)
            {
                Green_store[idx_e + st1] += -cimag(green_C[st1]) * kweight[kp];
                /*printf (" checkit  = %d %d %f \n", iene, idx + st1, Green_store[idx + st1]);*/
            }
        }
        idx_e += ntot;

    }                           /*  end for iene */



    /*===================================================================*/



    for (iene = 0; iene < EPoints; iene++)
    {

        root_pe = iene % pmo.npe_energy;
        idx = iene / pmo.npe_energy;


        for (st1 = 0; st1 < ntot; st1++)
        {
            density_matrix[st1] = Green_store[idx * ntot + st1];
        }


        idx = ntot;
        MPI_Bcast (density_matrix, idx, MPI_DOUBLE, root_pe,
                COMM_EN1);

        tri_to_row (density_matrix, work_matrix, ct.num_blocks, ct.block_dim);
        GetNewRho_on_c(states, rho, work_matrix);

        for (ix = 0; ix < get_FPX0_GRID(); ix++)
        {
		for (iy = 0; iy < get_FPY0_GRID(); iy++)
		{
			for (iz = 0; iz < get_FPZ0_GRID(); iz++)
			{
				idx = iz + iy * get_FPZ0_GRID() + ix * FPYZ;
				rho_energy[(ix + xoff)* get_FNY_GRID() * get_FNZ_GRID() + (iy + yoff)*  get_FNZ_GRID() + (iz + zoff)] += rho[idx];
			}
		}
        }
        MPI_Barrier(pct.img_comm);
    }



    iene = get_FNX_GRID() * get_FNY_GRID() * get_FNZ_GRID();
    global_sums (rho_energy, &iene, pct.grid_comm);
    if (pct.gridpe == 0)
    {
        double dx = get_celldm(0) / get_NX_GRID();
        double dy = get_celldm(0) * get_celldm(1) / get_NY_GRID();
        double dz = get_celldm(0) * get_celldm(2) / get_NZ_GRID();
        double B_A = 0.52917721;
        int count = 0;
        int level = 1; 
        char output[80];
	sprintf(output, "%s%d%s", "3Ddos_", number, ".cube");
        file = fopen (output, "w"); //create gaussian file to plot in PYMOL the 3D charge density for energy interval with high peak transmission
        fprintf( file, "Cubfile created from PWScf calculation\n" );
        fprintf( file, "Total SCF Density at Energy %12.9f \n", (EMIN+EMAX)/2 );
        fprintf( file, "1     0.000000    0.000000    0.000000 \n" );//hack the cube file by pretending there is only one atom in the gaussian file
        fprintf (file, "%d    %12.9f      0.000000    0.000000 \n", get_NX_GRID()/level, level*dx );//dx is the grid spacing in x in bohr
        fprintf (file, "%d    0.000000    %12.9f      0.000000 \n", get_NY_GRID()/level, level*dy );
        fprintf (file, "%d    0.000000    0.000000    %12.9f   \n", get_NZ_GRID()/level, level*dz );
        fprintf (file, "6     6.000000   10.000000   10.000000   10.000000 \n");//hack file by assigning just one carbon atom at some random position

        for (ix = 0; ix < get_NX_GRID()/level; ix++)
        {
            for (iy = 0; iy < get_NY_GRID()/level; iy++)
            {
                for (iz = 0; iz < get_NZ_GRID()/level; iz++)
                {
                    count ++;
                    fprintf ( file , " %18.6e",
                            rho_energy[ix * level * get_FG_RATIO() * get_FNY_GRID() * get_FNZ_GRID() + iy * level * get_FG_RATIO() *  get_FNZ_GRID() + iz * level * get_FG_RATIO()] );
                    if (count % 6 == 0)
                        fprintf (file, "\n");
                }
            }
        }


	fclose (file);
    }

    MPI_Barrier(pct.img_comm);
    fflush (NULL);

    /*===============================*/
//    my_free(tot);
//    my_free(tott);
//    my_free(g);
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
 //   my_free(ch0);
//    my_free(ch01);
//    my_free(ch10);
//    my_free(H10);
//    my_free(S10);



}
