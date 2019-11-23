/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "rmgthreads.h"
#include "vhartree.h"
#include "packfuncs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "Functional.h"
#include "Solvers.h"
#include "RmgParallelFft.h"

#include "blas.h"
#include "prototypes_tddft.h"
#include "RmgException.h"



template void RmgTddft<double> (double *, double *, double *,
          double *, double *, double *, double *, Kpoint<double> **);
template void RmgTddft<std::complex<double> > (double *, double *, double *,
          double *, double *, double *, double *, Kpoint<std::complex<double>> **);
template <typename OrbitalType> void RmgTddft (double * vxc, double * vh, double * vnuc, 
        double * rho, double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr)
{

    double *vtot, *vtot_psi, *vxc_psi=NULL;

    int dimx = Rmg_G->get_PX0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimy = Rmg_G->get_PY0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimz = Rmg_G->get_PZ0_GRID(Rmg_G->get_default_FG_RATIO());
    int FP0_BASIS = dimx * dimy * dimz;

    FILE *dfi = NULL;
    std::string filename;
    int n2,n22, numst, P0_BASIS,i, ione =1;
    int tot_steps = 0, pre_steps, tddft_steps;
    int Ieldyn = 1, iprint = 0;


    P0_BASIS =  Rmg_G->get_P0_BASIS(1);
    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    numst = ct.num_states; 
    n2 = numst * numst;
    n22 = 2* n2;

    double *Hmatrix = new double[n2];
    double *Hmatrix_old = new double[n2];
    double *Smatrix = new double[n2];
    double *Akick = new double[n2];
    double *Pn0 = new double[2*n2];
    double *Pn1 = new double[2*n2];
    double *vh_old = new double[FP0_BASIS];
    double *vxc_old = new double[FP0_BASIS];
    double *vh_corr_old = new double[FP0_BASIS];
    double *vh_corr = new double[FP0_BASIS];
//    double *vh_x = new double[FP0_BASIS];
//    double *vh_y = new double[FP0_BASIS];
//    double *vh_z = new double[FP0_BASIS];
    double *xpsi = new double[P0_BASIS * numst];
    double dipole_ele[3];


    RmgTimer *RT0 = new RmgTimer("2-TDDFT");

    // Loop over k-points
    if(ct.num_kpts != 1) 
    {
        rmg_printf(" \n  TDDFT does not support multiple k-points \n");

        fflush(NULL);
        throw RmgFatalException() << " TDDFT does not support multiple k-points in "<< __FILE__ << " at line " << __LINE__ << "\n";
        exit(0);
    }
    if(!ct.norm_conserving_pp)
    {
        rmg_printf(" \n  TDDFT support NCPP only \n");

        fflush(NULL);
        throw RmgFatalException() << " TDDFT support NCPP only in "<< __FILE__ << " at line " << __LINE__ << "\n";
        exit(0);
    }

    double efield[3];
    efield[0] = ct.x_field_0 * ct.e_field;
    efield[1] = ct.y_field_0 * ct.e_field;
    efield[2] = ct.z_field_0 * ct.e_field;
    if(pct.gridpe == 0)
    {
        filename = std::string(pct.image_path[pct.thisimg]) +"dipole.dat_"+ std::string(ct.basename);

        dfi = fopen(filename.c_str(), "w");

        fprintf(dfi, "\n  &&electric field:  %f  %f  %f ",efield[0], efield[1], efield[2]);

    }



    //    VhcorrDipoleInit(vh_x, vh_y, vh_z, rhoc);

    /* allocate memory for eigenvalue send array and receive array */

    vtot = new double[FP0_BASIS];
    vtot_psi = new double[P0_BASIS];
    if(ct.noncoll) vxc_psi = new double[4*P0_BASIS];
    
    double time_step =ct.tddft_time_step;

    if(ct.restart_tddft)
    {

        ReadData_rmgtddft(ct.outfile_tddft, vh, vxc, vh_corr, Pn0, Hmatrix, Smatrix,Smatrix, &pre_steps);
        dcopy(&n2, Hmatrix, &ione, Hmatrix_old, &ione);
        ReadData (ct.infile, vh, rho, vxc, Kptr);

    }
    else
    {
        for (int idx = 0; idx < FP0_BASIS; idx++) vtot[idx] = 0.0;
        init_efield(vtot);
        GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);

        HmatrixUpdate(Kptr[0], vtot_psi, (OrbitalType *)Akick);

        /* save old vhxc + vnuc */
        for (int idx = 0; idx < FP0_BASIS; idx++) {
            vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];
        }

        // Transfer vtot from the fine grid to the wavefunction grid
        GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);
        if(ct.noncoll)
        {
            for(int is = 0; is < 4; is++)
                GetVtotPsi (&vxc_psi[is*P0_BASIS], &vxc[is*FP0_BASIS], Rmg_G->default_FG_RATIO);
        }

        /*Generate the Dnm_I */
        get_ddd (vtot, vxc);

        HSmatrix (Kptr[0], vtot_psi, vxc_psi, (OrbitalType *)Hmatrix, (OrbitalType *)Smatrix);

        dcopy(&n2, Hmatrix, &ione, Hmatrix_old, &ione);

        if(pct.gridpe == 0 && ct.verbose)
        { 
            printf("\nHMa\n");
            for(i = 0; i < 10; i++) 
            {

                printf("\n");
                for(int j = 0; j < 10; j++) printf(" %8.1e",  Hmatrix[i*numst + j]);
            }

            printf("\nSMa\n");
            for(i = 0; i < 10; i++) 
            {

                printf("\n");
                for(int j = 0; j < 10; j++) printf(" %8.1e",  Smatrix[i*numst + j]);
            }
        }


        pre_steps = 0;
        for(i = 0; i < n2; i++) Hmatrix[i] += Akick[i]/time_step;

        for(i = 0; i < 2* n2; i++) Pn0[i] = 0.0;

        for(i = 0; i < ct.nel/2; i++) Pn0[i * numst + i] = 2.0;


        get_dipole(rho, dipole_ele);

        rmg_printf("\n  x dipolll  %f ", dipole_ele[0]);
        rmg_printf("\n  y dipolll  %f ", dipole_ele[1]);
        rmg_printf("\n  z dipolll  %f ", dipole_ele[2]);

 //         for(int i = 0; i < 10; i++) 
 //         { printf("Akick\n");
 //        for(int j = 0; j < 10; j++) printf(" %8.1e", i, Akick[i*numst + j]);
 //       }

    }



    for(tddft_steps = 0; tddft_steps < ct.tddft_steps; tddft_steps++)
    {

        tot_steps = pre_steps + tddft_steps;
        RmgTimer *RT2a = new RmgTimer("2-TDDFT: ELDYN");
        dscal(&n2, &time_step, Hmatrix, &ione);
        eldyn_(&numst, Smatrix, Hmatrix, Pn0, Pn1, &Ieldyn, &iprint);
        dcopy(&n22, Pn1, &ione, Pn0, &ione);
        delete(RT2a);
//            for(i = 0; i < 10; i++) 
//            { printf("Pn\n");
//           for(int j = 0; j < 10; j++) printf(" %8.1e", i, Pn1[i*numst + j]);
//          }


        RT2a = new RmgTimer("2-TDDFT: Rho");
        GetNewRho_rmgtddft((double *)Kptr[0]->orbital_storage, xpsi, rho, Pn1, numst);
        delete(RT2a);
        get_dipole(rho, dipole_ele);

        if(pct.gridpe == 0)fprintf(dfi, "\n  %f  %18.10f  %18.10f  %18.10f ",
                tot_steps*time_step, dipole_ele[0], dipole_ele[1], dipole_ele[2]);


        dcopy(&FP0_BASIS, vh_corr, &ione, vh_corr_old, &ione);
        dcopy(&FP0_BASIS, vh, &ione, vh_old, &ione);
        dcopy(&FP0_BASIS, vxc, &ione, vxc_old, &ione);

        //get_vxc(rho, rho_oppo, rhocore, vxc);
        double vtxc, etxc;
        RmgTimer *RT1 = new RmgTimer("2-TDDFT: exchange/correlation");
        Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
        F->v_xc(rho, rhocore, etxc, vtxc, vxc, ct.spin_flag );
        delete F;
        delete RT1;

        RT1 = new RmgTimer("2-TDDFT: Vh");
        VhDriver(rho, rhoc, vh, ct.vh_ext, 1.0-12);
        delete RT1;
        for (int idx = 0; idx < FP0_BASIS; idx++) {
            vtot[idx] = vxc[idx] + vh[idx]
                -vxc_old[idx] -vh_old[idx];
        }

        GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);
        RT2a = new RmgTimer("2-TDDFT: Hupdate");
        HmatrixUpdate(Kptr[0], vtot_psi, (OrbitalType *)Hmatrix);
        delete(RT2a);

        for(i = 0; i < n2; i++) Hmatrix[i] += Hmatrix_old[i];
        dcopy(&n2, Hmatrix, &ione, Hmatrix_old, &ione);

        if((tddft_steps +1) % ct.checkpoint == 0)
        {
            RT2a = new RmgTimer("2-TDDFT: Write");
            WriteData_rmgtddft(ct.outfile_tddft, vh, vxc, vh_corr, Pn0, Hmatrix, Smatrix, Smatrix, tot_steps);
            delete RT2a;
            fflush(NULL);
        }

    }

    if(pct.gridpe == 0) fclose(dfi);


    RmgTimer *RT2a = new RmgTimer("2-TDDFT: Write");
    WriteData_rmgtddft(ct.outfile_tddft, vh, vxc, vh_corr, Pn0, Hmatrix, Smatrix, Smatrix, tot_steps+1);
    delete RT2a;
    delete RT0;
}
