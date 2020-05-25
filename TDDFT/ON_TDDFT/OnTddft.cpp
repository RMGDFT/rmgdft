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
#include "LocalObject.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "Kbpsi.h"
#include "GpuAlloc.h"
#include "RmgException.h"
#include "blas_driver.h"
void eldyn_ort(int *desca, int Mdim, int Ndim, double *F,double *Po0,double *Po1,int *p_Ieldyn,  double *thrs,int*maxiter,  double *errmax,int
*niter , int *p_iprint, MPI_Comm comm) ;

void eldyn_nonort(int *p_N, double *S, double *F,double *Po0,double *Pn1,int *p_Ieldyn,  double *thrs,int*maxiter,  double *errmax,int *niter , int *p_iprint) ;



template void OnTddft<double>(double *, double *, double *,
          double *, double *, double *, double *, LocalObject<double> &, 
          LocalObject<double> &, LocalObject<double> &); 
//template void OnTddft<std::complex<double> >(double *, double *, double *,
//        double *, double *, double *, double *, LocalObject<std::complex<double>> &,
//        LocalObject<std::complex<double>> &,LocalObject<std::complex<double>> &);
template <typename OrbitalType> void OnTddft (double * vxc, double * vh, double * vnuc, 
        double * rho, double * rho_oppo, double * rhocore, double * rhoc, LocalObject<OrbitalType> &Phi,
        LocalObject<OrbitalType> &H_Phi, LocalObject<OrbitalType> &LP)
{

    double *vtot, *vtot_psi;

    int dimx = Rmg_G->get_PX0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimy = Rmg_G->get_PY0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimz = Rmg_G->get_PZ0_GRID(Rmg_G->get_default_FG_RATIO());
    int FP0_BASIS = dimx * dimy * dimz;

    FILE *dfi = NULL;
    std::string filename;
    int n2,n22, numst, P0_BASIS,i, ione =1;
    int tot_steps=0, pre_steps, tddft_steps;
    int Ieldyn = 1, iprint = 0;
    double one = 1.0, zero = 0.0;


    P0_BASIS =  Rmg_G->get_P0_BASIS(1);
    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    numst = Phi.num_tot;
    int scalapack_groups = 1;
#if GPU_ENABLED
    scalapack_groups = pct.grid_npes;
#endif
    int last = 1;
    numst = ct.num_states;
    Scalapack *Sp = new Scalapack(scalapack_groups, pct.thisimg, ct.images_per_node, numst,
            ct.scalapack_block_factor, last, pct.grid_comm);
    int Mdim = Sp->GetDistMdim();
    int Ndim = Sp->GetDistNdim();
    n2 = Sp->GetDistMdim() * Sp->GetDistNdim();

    if(!Sp->Participates()) n2 = 1;
    int *desca = Sp->GetDistDesca();

    n22 = 2* n2;

    double *Hmatrix_glob = (double *)GpuMallocManaged((size_t)numst * (size_t)numst*sizeof(double));
    double *Smatrix_glob = (double *)GpuMallocManaged((size_t)numst * (size_t)numst*sizeof(double));
    double *Smatrix     = (double *)GpuMallocManaged((size_t)n2*sizeof(double));
    double *Hmatrix     = (double *)GpuMallocManaged((size_t)n2*sizeof(double));
    double *Hmatrix_old = (double *)GpuMallocManaged((size_t)n2*sizeof(double));
    double *Akick       = (double *)GpuMallocManaged((size_t)n2*sizeof(double));
    double *Pn0         = (double *)GpuMallocManaged((size_t)n2*sizeof(double)*2);
    double *Pn1         = (double *)GpuMallocManaged((size_t)n2*sizeof(double)*2);
    double *vh_old      = new double[FP0_BASIS];
    double *vxc_old     = new double[FP0_BASIS];
    double *vh_corr_old = new double[FP0_BASIS];
    double *vh_corr     = new double[FP0_BASIS];
    // Jacek:
    //double *dHmatrix    = new double[n2];   // storage for  H1 -H1_old
    double *Hmatrix_m1  = (double *)GpuMallocManaged((size_t)n2*sizeof(double));
    double *Hmatrix_0   = (double *)GpuMallocManaged((size_t)n2*sizeof(double));
    double *Hmatrix_1   = (double *)GpuMallocManaged((size_t)n2*sizeof(double));
    double *Hmatrix_dt  = (double *)GpuMallocManaged((size_t)n2*sizeof(double));
    double    err        ;
    int       ij_err     ;
    double   thrs_dHmat =1.0e-5 ;

    if(pct.gridpe == 0) {
        printf("\n Number of states used for TDDFT: Nbasis =  %d \n",numst);
        printf(" Propagator used :  Ieldyn = %d  \\1=BCH, 2=Diagonalizer\\ \n",Ieldyn) ;
    }

    double *Cmatrix     = (double *)GpuMallocManaged((size_t)n2*sizeof(double));
    double *Xmatrix     = (double *)GpuMallocManaged((size_t)n2*sizeof(double));

    //    double *vh_x = new double[FP0_BASIS];
    //    double *vh_y = new double[FP0_BASIS];
    //    double *vh_z = new double[FP0_BASIS];

    int num_orb = Phi.num_thispe;
    double *rho_matrix_local = (double *)GpuMallocManaged(num_orb * num_orb * sizeof(double));
    double *Hij_local = (double *)GpuMallocManaged(num_orb * num_orb * sizeof(double));
    double *Sij_local = (double *)GpuMallocManaged(num_orb * num_orb * sizeof(double));
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
        filename = std::string(pct.image_path[pct.thisimg])+ "dipole.dat_"+std::string(ct.basename);

        dfi = fopen(filename.c_str(), "w");

        fprintf(dfi, "\n  &&electric field:  %f  %f  %f ",efield[0], efield[1], efield[2]);

    }



    //    VhcorrDipoleInit(vh_x, vh_y, vh_z, rhoc);

    /* allocate memory for eigenvalue send array and receive array */

    vtot = new double[FP0_BASIS];
    vtot_psi = new double[P0_BASIS];
    double time_step =ct.tddft_time_step;

    if(ct.restart_tddft)
    {

        ReadData_rmgtddft(ct.outfile_tddft, vh, vxc, vh_corr, Pn0, Hmatrix, Smatrix, 
                Cmatrix, Hmatrix_m1, Hmatrix_0, &pre_steps);
        dcopy(&n2, Hmatrix, &ione, Hmatrix_old, &ione);
        Phi.ReadOrbitalsFromSingleFiles(std::string(ct.infile), *Rmg_G);
    }
    else
    {

        /* save old vhxc + vnuc */
        for (int idx = 0; idx < FP0_BASIS; idx++) {
            vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];
        }

        // Transfer vtot from the fine grid to the wavefunction grid
        GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);

        /*Generate the Dnm_I */
        get_ddd (vtot, vxc, true);

        LO_x_LO(LP, Phi, Kbpsi_mat_local, *Rmg_G);
        mat_local_to_glob(Kbpsi_mat_local, Kbpsi_mat, LP, Phi, 0, LP.num_tot, 0, Phi.num_tot, true);
        ApplyHphi(Phi, H_Phi, vtot_psi);

        LO_x_LO(Phi, Phi, Sij_local, *Rmg_G);

        LO_x_LO(Phi, H_Phi, Hij_local, *Rmg_G);

        mat_local_to_glob(Sij_local, Smatrix_glob, Phi, Phi, 0, Phi.num_tot, 0, Phi.num_tot, false);
        mat_local_to_glob(Hij_local, Hmatrix_glob, Phi, Phi, 0, Phi.num_tot, 0, Phi.num_tot, false);



        GetHvnlij_proj(Hmatrix_glob, Smatrix_glob, Kbpsi_mat, Kbpsi_mat,
                Phi.num_tot, Phi.num_tot, LP.num_tot, true);

        int idx = Phi.num_tot * Phi.num_tot;
        MPI_Allreduce(MPI_IN_PLACE, Hmatrix_glob, idx, MPI_DOUBLE, MPI_SUM, Phi.comm);
        MPI_Allreduce(MPI_IN_PLACE, Smatrix_glob, idx, MPI_DOUBLE, MPI_SUM, Phi.comm);
        Sp->CopySquareMatrixToDistArray(Hmatrix_glob, Hmatrix, numst, desca);
        Sp->CopySquareMatrixToDistArray(Smatrix_glob, Smatrix, numst, desca);

//        DiagScalapack(states, ct.num_states, Hmatrix, Smatrix);

#if GPU_ENABLED
        mat_dist_to_global(zz_dis, pct.desca, Cmatrix);
#else
        dcopy(&n2, zz_dis, &ione, Cmatrix, &ione);
#endif
        dgemm_driver ("T", "N", numst, numst, numst, one, Cmatrix, ione, ione, desca,
            Hmatrix, ione, ione, desca, zero, Akick, ione, ione, desca);
        dgemm_driver ("N", "N", numst, numst, numst, one, Akick, ione, ione, desca,
            Cmatrix, ione, ione, desca, zero, Hmatrix_old, ione, ione, desca);
        dgemm_driver ("T", "N", numst, numst, numst, one, Cmatrix, ione, ione, desca,
            Smatrix, ione, ione, desca, zero, Akick, ione, ione, desca);
        dgemm_driver ("N", "N", numst, numst, numst, one, Akick, ione, ione, desca,
            Cmatrix, ione, ione, desca, zero, Smatrix, ione, ione, desca);
        my_sync_device();

        if(pct.gridpe == 0 && ct.verbose)
        { 
            printf("\nHMa\n");
            for(i = 0; i < 10; i++) 
            {

                printf("\n");
                for(int j = 0; j < 10; j++) printf(" %8.1e",  Hmatrix_old[i*numst + j]);
            }

            printf("\nSMa\n");
            for(i = 0; i < 10; i++) 
            {

                printf("\n");
                for(int j = 0; j < 10; j++) printf(" %8.1e",  Smatrix[i*numst + j]);
            }
        }


        pre_steps = 0;

        for (int idx = 0; idx < FP0_BASIS; idx++) vtot[idx] = 0.0;
        init_efield(vtot);
        GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);

        HmatrixUpdate_on(Phi, H_Phi, vtot_psi, Hij_local);
        mat_local_to_glob(Hij_local, Hmatrix_glob, Phi, Phi, 0, Phi.num_tot, 0, Phi.num_tot, true);
        Sp->CopySquareMatrixToDistArray(Hmatrix_glob, Akick, numst, desca);

        double alpha = 1.0/time_step;
        daxpy_driver ( n2 ,  alpha, Akick, ione , Hmatrix ,  ione) ;

        for(i = 0; i < 2* n2; i++) Pn0[i] = 0.0;

        for(i = 0; i < numst * numst; i++) Hmatrix_glob[i] = 0.0;
        for(i = 0; i < ct.nel/2; i++)  Hmatrix_glob[i * numst + i] = 2.0;
        Sp->CopySquareMatrixToDistArray(Hmatrix_glob, Pn0, numst, desca);


        get_dipole(rho, dipole_ele);

        rmg_printf("\n  x dipolll  %f ", dipole_ele[0]);
        rmg_printf("\n  y dipolll  %f ", dipole_ele[1]);
        rmg_printf("\n  z dipolll  %f ", dipole_ele[2]);

        //         for(int i = 0; i < 10; i++) 
        //         { printf("Akick\n");
        //        for(int j = 0; j < 10; j++) printf(" %8.1e", i, Akick[i*numst + j]);
        //       }
        dgemm_driver ("T", "N", numst, numst, numst, one, Cmatrix, ione, ione, desca,
                Hmatrix, ione, ione, desca, zero, Akick, ione, ione, desca);
        dgemm_driver ("N", "N", numst, numst, numst, one, Akick, ione, ione, desca,
                Cmatrix, ione, ione, desca, zero, Hmatrix, ione, ione, desca);

        if(pct.gridpe == 0 && ct.verbose)
        { 
            printf("\nHMc\n");
            for(i = 0; i < 10; i++) 
            {

                printf("\n");
                for(int j = 0; j < 10; j++) printf(" %8.1e",  Hmatrix[i*MXLLDA + j]);
            }
        }

        dcopy_driver(n2, Hmatrix, ione, Hmatrix_m1, ione);
        dcopy_driver(n2, Hmatrix, ione, Hmatrix_0 , ione);
        my_sync_device();

    }




    for(tddft_steps = 0; tddft_steps < ct.tddft_steps; tddft_steps++)
    {

        tot_steps = pre_steps + tddft_steps;

        extrapolate_Hmatrix  (Hmatrix_m1,  Hmatrix_0, Hmatrix_1  , n2) ; //   (*Hm1, double *H0, double *H1,  int *ldim)

        my_sync_device();
        //  SCF loop 
        int  Max_iter_scf = 10 ; int  iter_scf =0 ;
        err =1.0e0   ;  thrs_dHmat  = 1e-7  ;

        double  thrs_bch =1.0e-10; 
        int     maxiter_bch  =100;
        double  errmax_bch ;
        int     niter_bch ;

        RmgTimer *RT2a ;    // timer type  declaration

        //-----   SCF loop  starts here: 
        while (err > thrs_dHmat &&  iter_scf <  Max_iter_scf)  {

            //RmgTimer *RT2a = new RmgTimer("2-TDDFT: ELDYN");
            RT2a = new RmgTimer("2-TDDFT: ELDYN");
            magnus (Hmatrix_0,    Hmatrix_1 , time_step, Hmatrix_dt , n2) ; 
            /* --- fortran version:  --*/
            // eldyn_(&numst, Smatrix, Hmatrix_dt, Pn0, Pn1, &Ieldyn, &iprint);
            /* --- C++  version:  --*/
            my_sync_device();
            eldyn_ort(desca, Mdim, Ndim,  Hmatrix_dt,Pn0,Pn1,&Ieldyn, &thrs_bch,&maxiter_bch,  &errmax_bch,&niter_bch ,  &iprint, Sp->GetComm()) ;

            delete(RT2a);

            //if(pct.gridpe == 0 && ct.verbose)
            //{ 
            //    printf("\nPn1\n");
            //    for(i = 0; i < 10; i++) 
            //    {
//
 //                   printf("\n");
  //                  for(int j = 0; j < 10; j++) printf(" %8.1e",  Pn1[i*MXLLDA + j]);
   //             }
    //        }

            RT2a = new RmgTimer("2-TDDFT: transform of H to Cij *H * Cij");
            dgemm_driver ("N", "N", numst, numst, numst, one, Cmatrix, ione, ione, desca,
                    Pn1, ione, ione, desca, zero, Akick, ione, ione, desca);
            dgemm_driver ("N", "T", numst, numst, numst, one, Akick, ione, ione, desca,
                    Cmatrix, ione, ione, desca, zero, Xmatrix, ione, ione, desca);

            if( scalapack_groups != pct.grid_npes)
            {
                if(Sp->Participates()) 
                {
                    Sp->CopyDistArrayToSquareMatrix(Hmatrix_glob, Xmatrix, numst, desca);
                    if( scalapack_groups != pct.grid_npes)
                        Sp->ScalapackBlockAllreduce(Hmatrix_glob, (size_t)numst * (size_t)numst);
                }

                Sp->BcastRoot(Hmatrix_glob, numst * numst, MPI_DOUBLE);
            }
            else
            {
                dcopy_driver(n2, Xmatrix, ione, Hmatrix_glob, ione);
            }


            my_sync_device();


            //      printf("\n PPP \n");
            //      for(i = 0; i < 10; i++) 
            //      { printf("\n");
            //          for(int j = 0; j < 10; j++) printf(" %8.2e", Pn1[i*numst + j]);
            //      }
            //      printf("\n\n");


            RT2a = new RmgTimer("2-TDDFT: mat_glob_to_local");

            mat_global_to_local(Phi, H_Phi, Hmatrix_glob, rho_matrix_local);
            delete(RT2a);
            RT2a = new RmgTimer("2-TDDFT: Rho");
            GetNewRho_proj(Phi, H_Phi, rho, rho_matrix_local);


            double tcharge = 0.0;
            for (i = 0; i < get_FP0_BASIS(); i++)
                tcharge += rho[i];
            ct.tcharge = real_sum_all(tcharge, pct.grid_comm);
            ct.tcharge = real_sum_all(ct.tcharge, pct.spin_comm);


            ct.tcharge *= get_vel_f();

            double t2 = ct.nel / ct.tcharge;
            dscal(&FP0_BASIS, &t2, rho, &ione);


            if(fabs(t2 -1.0) > 1.0e-11 && pct.gridpe == 0)
                rmg_printf("\n Warning: total charge Normalization constant = %e  \n", t2-1.0);

            delete(RT2a);


            dcopy(&FP0_BASIS, vh_corr, &ione, vh_corr_old, &ione);
            dcopy(&FP0_BASIS, vh, &ione, vh_old, &ione);
            dcopy(&FP0_BASIS, vxc, &ione, vxc_old, &ione);

            //get_vxc(rho, rho_oppo, rhocore, vxc);
            double vtxc, etxc;
            RmgTimer *RT1 = new RmgTimer("2-TDDFT: exchange/correlation");
            Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
            F->v_xc(rho, rhocore, etxc, vtxc, vxc, ct.nspin );
            delete F;
            delete RT1;

            RT1 = new RmgTimer("2-TDDFT: Vh");
            VhDriver(rho, rhoc, vh, ct.vh_ext, 1.0-12);
            delete RT1;
            double tem = 0.0;
            for (int idx = 0; idx < FP0_BASIS; idx++) {
                vtot[idx] = vxc[idx] + vh[idx]
                    -vxc_old[idx] -vh_old[idx];
                tem += vtot[idx] * vtot[idx];
            }

            MPI_Allreduce(MPI_IN_PLACE, &tem, 1, MPI_DOUBLE, MPI_SUM, Phi.comm);
            if(pct.gridpe == 0) printf(" pote diff  %e", tem);
            RT1 = new RmgTimer("2-TDDFT: Hupdate");
            GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);
            HmatrixUpdate_on(Phi, H_Phi, vtot_psi, Hij_local);
            mat_local_to_glob(Hij_local, Hmatrix_glob, Phi, Phi, 0, Phi.num_tot, 0, Phi.num_tot, true);

            Sp->CopySquareMatrixToDistArray(Hmatrix_glob, Hmatrix, numst, desca);
            dgemm_driver ("T", "N", numst, numst, numst, one, Cmatrix, ione, ione, desca,
                    Hmatrix, ione, ione, desca, zero, Akick, ione, ione, desca);
            dgemm_driver ("N", "N", numst, numst, numst, one, Akick, ione, ione, desca,
                    Cmatrix, ione, ione, desca, zero, Hmatrix, ione, ione, desca);


            delete RT1;

            my_sync_device();
            double one = 1.0;
            daxpy_driver ( n2 ,  one, Hmatrix_old, ione , Hmatrix ,  ione) ;
            dcopy_driver(n2, Hmatrix, ione, Hmatrix_old, ione);         // saves Hmatrix to Hmatrix_old   

            //////////  < ---  end of Hamiltonian update

            // check error and update Hmatrix_1:
            my_sync_device();
            tst_conv_matrix (&err, &ij_err ,  Hmatrix,  Hmatrix_1 ,  n2, Sp->GetComm()) ;  //  check error  how close  H and H_old are
            dcopy_driver(n2, Hmatrix  , ione, Hmatrix_1, ione);

            if(pct.gridpe == 0) { printf("step: %5d  iteration: %d  thrs= %12.5e err=  %12.5e at element: %5d \n", 
                    tddft_steps, iter_scf,    thrs_dHmat,  err,         ij_err); } 
            //err= -1.0e0 ;  
            iter_scf ++ ;
        } //---- end of  SCF/while loop 

        /*  done with propagation,  save Pn1 ->  Pn0 */
        dcopy_driver(n22, Pn1, ione, Pn0, ione);

        //  extract dipole from rho(Pn1)
        get_dipole(rho, dipole_ele);
        // save current  H0, H1 for the  next step extrapolatiion
        dcopy_driver  (n2, Hmatrix_0, ione, Hmatrix_m1 , ione);         
        //dcopy(&n2, Hmatrix  , &ione, Hmatrix_1  , &ione);         // this update is already done right after scf loop 

        dcopy_driver  (n2, Hmatrix_1, ione, Hmatrix_0  , ione);        

        if(pct.gridpe == 0)fprintf(dfi, "\n  %f  %18.10f  %18.10f  %18.10f ",
                tot_steps*time_step, dipole_ele[0], dipole_ele[1], dipole_ele[2]);


        if((tddft_steps +1) % ct.checkpoint == 0)
        {
            RmgTimer *RT1 = new RmgTimer("2-TDDFT: WriteData");
            my_sync_device();
            WriteData_rmgtddft(ct.outfile_tddft, vh, vxc, vh_corr, Pn0, Hmatrix, Smatrix,
                    Cmatrix, Hmatrix_m1, Hmatrix_0, tot_steps);
            delete RT1;
            fflush(NULL);
        }

    }

    if(pct.gridpe == 0) fclose(dfi);


    my_sync_device();
    WriteData_rmgtddft(ct.outfile_tddft, vh, vxc, vh_corr, Pn0, Hmatrix, Smatrix, 
            Cmatrix, Hmatrix_m1, Hmatrix_0, tot_steps+1);
    delete RT0;
}
