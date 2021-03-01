#include "negf_prototypes.h"
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
#include <complex>

//#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"
#include "prototypes_on.h"
#include "GpuAlloc.h"
#include "prototypes_negf.h"
#include "transition.h"
#include "my_mpi.h"

void get_dos (STATE * states)
{
    int iprobe, idx_e;
    int iene, st1;
    std::complex<double> *sigma;
    double de, emin, emax;

    std::complex<double> *work, *green_C;
    std::complex<double> ene;
    std::complex<double> I(0.0, 1.0);
    int kp, i, *sigma_idx, idx_C;
    FILE *file;

    int ntot, ndim;
    int  xoff, yoff, zoff;
    int root_pe, idx, ix, iy, iz;

    int E_POINTS, nkp[3];
    double E_imag, KT;
    int FPYZ = get_FPY0_GRID() * get_FPZ0_GRID();
    int nx1, nx2, ny1, ny2, nz1, nz2; 
    double *kvecx, *kvecy, *kvecz, *kweight;


    read_cond_input (&emin, &emax, &E_POINTS, &E_imag, &KT, nkp);
    de = (emax - emin) / (E_POINTS - 1);

    //  set number of kpoint in x direction be 1, 
    nkp[0] = 1;

    if(cei.num_probe > 2 ) nkp[1] = 1;
    int nkp_tot = nkp[0] * nkp[1] * nkp[2];
    rmg_printf("\n nkp  %d %d %d", nkp[0], nkp[1], nkp[2]);

    if (nkp_tot == 0 ) rmg_error_handler (__FILE__, __LINE__, "wrong number of kpoints in cond.in");

    kvecx = new double[nkp_tot];
    kvecy = new double[nkp_tot];
    kvecz = new double[nkp_tot];
    kweight = new double[nkp_tot];

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
    lcr[0].Htri = (double *)RmgMallocHost(pmo.ntot * sizeof(double));
    lcr[0].Stri = (double *)RmgMallocHost(pmo.ntot * sizeof(double));

    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx = pmo.mxllda_lead[iprobe-1] * pmo.mxlocc_lead[iprobe-1];
        lcr[iprobe].H00 = (double *)RmgMallocHost(idx * sizeof(double));
        lcr[iprobe].S00 = (double *)RmgMallocHost(idx * sizeof(double));
        lcr[iprobe].H01 = (double *)RmgMallocHost(idx * sizeof(double));
        lcr[iprobe].S01 = (double *)RmgMallocHost(idx * sizeof(double));
    }

    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        i = cei.probe_in_block[iprobe - 1];
        idx = pmo.mxllda_cond[i] * pmo.mxlocc_lead[iprobe-1];
        lcr[iprobe].HCL = (double *)RmgMallocHost(idx * sizeof(double));
        lcr[iprobe].SCL = (double *)RmgMallocHost(idx * sizeof(double));
    }


    read_matrix_pp();


    /*======================== Allocate memory for sigma =================*/

    idx = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        idx = std::max(idx, pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]);
    }
    sigma = (std::complex<double> *)RmgMallocHost(idx * sizeof(std::complex<double>));

    sigma_idx = new int[cei.num_probe];

    idx = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        sigma_idx[iprobe - 1] = idx;
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        idx += pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C];
    }

    sigma_all = (std::complex<double> *)RmgMallocHost(idx * sizeof(std::complex<double>));


    /*============== Allocate memory for tot, tott, g ====================*/

    idx = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
        idx = std::max(idx, pmo.mxllda_cond[idx_C] * pmo.mxlocc_lead[iprobe-1]);
    }

    work = (std::complex<double> *)RmgMallocHost(12*idx * sizeof(std::complex<double>));
    green_C = (std::complex<double> *)RmgMallocHost(pmo.ntot_low * sizeof(std::complex<double>));
    st1 = ( E_POINTS + pmo.npe_energy-1)/pmo.npe_energy;
    double *Green_store = (double *)RmgMallocHost(st1 * pmo.ntot * sizeof(double));
    for(int idx = 0; idx < st1 * pmo.ntot; idx++) Green_store[idx] = 0.0;


    double *density_matrix =(double *)RmgMallocHost(pmo.ntot * sizeof(double));

    /*===================================================================*/

    nx1 = cei.dos_window_start[0] * get_FG_RATIO();
    nx2 = cei.dos_window_end[0] * get_FG_RATIO();
    ny1 = cei.dos_window_start[1] * get_FG_RATIO();
    ny2 = cei.dos_window_end[1] * get_FG_RATIO();
    nz1 = cei.dos_window_start[2] * get_FG_RATIO();
    nz2 = cei.dos_window_end[2] * get_FG_RATIO();


    double *rho_energy =(double *)RmgMallocHost(E_POINTS * get_FNX_GRID() * sizeof(double));
    double *rho_energy2 =(double *)RmgMallocHost(E_POINTS * get_FNY_GRID() * sizeof(double));
    for(int idx = 0; idx< E_POINTS * get_FNX_GRID(); idx++) rho_energy[idx] = 0.0;
    for(int idx = 0; idx< E_POINTS * get_FNY_GRID(); idx++) rho_energy2[idx] = 0.0;



    xoff = get_FPX_OFFSET();
    yoff = get_FPY_OFFSET();
    zoff = get_FPZ_OFFSET();


    /*===================================================================*/


    idx_e = 0;

    /*for (iene = pct.gridpe; iene < E_POINTS; iene += pct.grid_npes)*/
    for (iene = pmo.myblacs; iene < E_POINTS; iene += pmo.npe_energy)
    {
        ene = emin + iene * de + I * E_imag;


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

                idx_C = cei.probe_in_block[iprobe - 1];  /* block index */
                for (i = 0; i < pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C]; i++)
                {
                    sigma_all[sigma_idx[iprobe - 1] + i] = sigma[i];
                }

            }                       /*  end for iprobe */


            Sgreen_c_p (lcr[0].Htri, lcr[0].Stri, sigma_all, sigma_idx, 
                    ene, green_C); 


            for (st1 = 0; st1 < ntot; st1++)
            {
                Green_store[idx_e + st1] += -std::imag(green_C[st1]) * kweight[kp];
                /*printf (" checkit  = %d %d %f \n", iene, idx + st1, Green_store[idx + st1]);*/
            }

        }

        idx_e += ntot;
    }                           /*  end for iene */




    /*===================================================================*/

    size_t size = LocalOrbital->num_thispe * LocalOrbital->num_thispe * sizeof(double);
    double *rho_matrix_local = (double *)RmgMallocHost(size);

    /*for (iene = pmo.myblacs; iene < E_POINTS; iene += pmo.npe_energy)*/
    for (iene = 0; iene < E_POINTS; iene++)
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

        RmgTimer *RT4 = new RmgTimer("3-SCF: rho");
        //    get_new_rho_soft (states, rho);
        tri_to_local (density_matrix, rho_matrix_local, *LocalOrbital);
        GetNewRho_proj(*LocalOrbital, *H_LocalOrbital, rho, rho_matrix_local);
        //    get_new_rho_local (states_distribute, rho);
        delete(RT4);


        for (ix = 0; ix < get_FPX0_GRID(); ix++)
        {
            for (iy = 0; iy < get_FPY0_GRID(); iy++)
            {
                if((iy + yoff) >= ny1 && (iy + yoff) <= ny2)
                {
                    for (iz = 0; iz < get_FPZ0_GRID(); iz++)
                    {
                        if((iz + zoff) >= nz1 && (iz + zoff) <= nz2)
                        {
                            idx = iz + iy * get_FPZ0_GRID() + ix * FPYZ;
                            rho_energy[iene * get_FNX_GRID() + ix + xoff] += rho[idx];
                        }
                    }
                }
            }
        }


        if (cei.num_probe > 2)
        {
            for (iy = 0; iy < get_FPY0_GRID(); iy++)
            {
                for (ix = 0; ix < get_FPX0_GRID(); ix++)
                {
                    if((ix + xoff) >= nx1 && (ix + xoff) <= nx2)
                    {
                        for (iz = 0; iz < get_FPZ0_GRID(); iz++)
                        {
                            if((iz + zoff) >= nz1 && (iz + zoff) <= nz2)
                            {
                                idx = iz + iy * get_FPZ0_GRID() + ix * FPYZ;
                                rho_energy2[iene * get_FNY_GRID() + iy + yoff] += rho[idx];
                            }
                        }
                    }
                }
            }
        }

        MPI_Barrier(pct.img_comm);
    }



    iene = E_POINTS * get_FNX_GRID();
    global_sums (rho_energy, &iene, pct.grid_comm);
    if (pct.gridpe == 0)
    {
        double dx = get_celldm(0) / get_NX_GRID();
        double x0 = 0.5 * get_celldm(0);

        file = fopen ("dos.dat", "w");
        fprintf (file, "#     x[a0]      E[eV]          dos\n\n");
        for (iene = 0; iene < E_POINTS; iene++)
        {

            for (ix = 0; ix < get_NX_GRID(); ix++)
            {

                fprintf (file, " %10.6f %10.6f %12.6e\n",
                        ix * dx - x0, emin+iene*de, rho_energy[iene * get_FNX_GRID() + ix * get_FG_RATIO()]);
            }
            fprintf (file, "\n");
        }

        fclose (file);
    }

    if (cei.num_probe > 2)
    {

        iene = E_POINTS * get_FNY_GRID();
        global_sums (rho_energy2, &iene, pct.grid_comm);
        if (pct.gridpe == 0)
        {
            double y = get_celldm(1) * get_celldm(0);
            double dy = y / get_NY_GRID();
            double y0 = 0.5 * y;

            file = fopen ("dos2.dat", "w");
            fprintf (file, "#     y[b0]      E[eV]          dos\n\n");
            for (iene = 0; iene < E_POINTS; iene++)
            {

                for (iy = 0; iy < get_NY_GRID(); iy++)
                {

                    fprintf (file, " %10.6f %10.6f %12.6e\n",
                            iy * dy - y0, emin+iene*de, rho_energy2[iene * get_FNY_GRID() + iy * get_FG_RATIO()]);
                }
                fprintf (file, "\n");
            }

            fclose (file);
        }
    }


    MPI_Barrier(pct.img_comm);
    fflush (NULL);

    /*===============================*/
    delete [] sigma_idx;
    RmgFreeHost(sigma_all);
    RmgFreeHost(green_C);

    RmgFreeHost(lcr[0].Htri);
    RmgFreeHost(lcr[0].Stri);
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        RmgFreeHost(lcr[iprobe].H00);
        RmgFreeHost(lcr[iprobe].S00);
        RmgFreeHost(lcr[iprobe].H01);
        RmgFreeHost(lcr[iprobe].S01);
    }

    /*===============================*/

    RmgFreeHost(rho_energy); 
    RmgFreeHost(rho_energy2); 
    RmgFreeHost(Green_store);
    RmgFreeHost(rho_matrix_local);


}

