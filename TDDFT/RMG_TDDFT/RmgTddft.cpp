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
#include "blas_driver.h"
#include "GpuAlloc.h"


void  init_point_charge_pot(double *vtot_psi, int density);
void eldyn_ort(int *desca, int Mdim, int Ndim, double *F,double *Po0,double *Po1,int *p_Ieldyn,  double *thrs,int*maxiter,  double *errmax,int
*niter , int *p_iprint, MPI_Comm comm) ;
void eldyn_nonort(int *p_N, double *S, double *F,double *Po0,double *Pn1,int *p_Ieldyn,  double *thrs,int*maxiter,  double *errmax,int *niter , int *p_iprint) ;
  

void print_matrix_d(double *matrix,  int  *nblock, int *ldim ){
    /*    
          Print the block [0:nblock-1,0:nblock-1]  of an  array of size
          ldim*ldim.  Array (matrix) is stored as a vector.

     */ 

    int  N  = *nblock ;   // block size to be printed 
    int  LD = *ldim   ;   // leaading dimension of matrix

    int i,j, ij   ;

    if (pct.gridpe ==0 ) {
        for     ( i = 0; i< N  ; i++){   // row
            //printf("\n   %d  : ", i);
            printf("\n"); 
            for  ( j = 0; j< N  ; j++){   // column
                ij = i*LD  + j  ;
                printf(" %14.6e", matrix[ij]) ;
            }
        }
        printf("\n");
    }
}

void write_rho_x(double * rho, char *ab)
{

    int ix, iy, iz, poff;
    double t1;
    double *zvec;

    zvec = new double[get_FNX_GRID()]();
    /* Get this processors offset */
    poff = get_FPX_OFFSET();



    /* Zero out result vector */
    for (ix = 0; ix < get_FNX_GRID(); ix++)
        zvec[ix] = ZERO;


    /* Loop over this processor */
    for (ix = 0; ix < get_FPX0_GRID(); ix++)
    {
        t1 = 0.0;
        for (iy = 0; iy < get_FPY0_GRID(); iy++)
            for (iz = 0; iz < get_FPZ0_GRID(); iz++)

                t1 += rho[ix * get_FPY0_GRID() * get_FPZ0_GRID() + iy * get_FPZ0_GRID() + iz];


        zvec[ix + poff] = t1;

    }                           /* end for */


    /* Now sum over all processors */
    ix = get_FNX_GRID();
    global_sums(zvec, &ix, pct.grid_comm);

    if (pct.gridpe == 0)
    {
        printf("\n\n Planar average of the electrostatic density\n");
        for (ix = 0; ix < get_FNX_GRID(); ix++)
        {
            t1 = ix * get_hxxgrid();
            printf(" %d %f %s\n", ix, zvec[ix] / get_FNY_GRID() / get_FNZ_GRID(), ab);
        }
        printf(" & %s\n", ab);
        fflush(NULL);
    }

    delete [] zvec;
}                               /* end get_avgd */

/******/
///////////////////////////////////////

void print_matrix_z(double *matrix,  int  *nblock, int *ldim ){
    /*    
          Print the block [0:nblock-1,0:nblock-1]  of a complex  array of size
          ldim*ldim.  Array (matrix) is stored as a ldim*ldim vector or real  
          values followed by ldim*ldim vector of imaginary values

     */ 

    int  N  = *nblock ;   // block size to be printed 
    int  LD = *ldim   ;   // leaading dimension of matrix

    int i,j, ij   ;
    int i0_re = 0 ; 
    int i0_im = i0_re + LD*LD  ; 

    if (pct.gridpe ==0 ) {

        // print  real part  
        printf("***  real:\n");
        for     ( i = 0; i< N  ; i++){         // row
            printf("\n");                      // printf("   %d  : ", i);
            for  ( j = 0; j< N  ; j++){        // column
                ij = i*LD  + j  ;
                printf(" %14.6e", matrix[i0_re+ ij]) ;
            }
        }
        printf("\n") ;   

        // print  imaginary  part  
        printf("***  imag:\n ") ;   
        for     ( i = 0; i< N  ; i++){         // row
            printf("\n");                       // printf("   %d  : ", i);
            for  ( j = 0; j< N  ; j++){         // column
                ij = i*LD  + j  ;
                printf(" %14.6e", matrix[ i0_im+ij ]) ;
            }
        }
        printf("\n");
    }
}


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
    int Ieldyn = 1;    // BCH  
    //int Ieldyn = 2;    // Diagev
    int iprint = 0;


    P0_BASIS =  Rmg_G->get_P0_BASIS(1);
    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    int scalapack_groups = 1;
    switch(ct.subdiag_driver) {

        case SUBDIAG_LAPACK:
            scalapack_groups = pct.grid_npes;
            break;
        case SUBDIAG_SCALAPACK:
            scalapack_groups = 1;
            if(ct.scalapack_block_factor >= ct.num_states)
                scalapack_groups = pct.grid_npes;
            break;
        case SUBDIAG_CUSOLVER:
            scalapack_groups = pct.grid_npes;
            break;
        default:
            rmg_error_handler(__FILE__, __LINE__, "Invalid subdiag_driver type in TDDFT");

    } // end switch

#if CUDA_ENABLED
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

    double *matrix_glob = (double *)GpuMallocManaged((size_t)numst * (size_t)numst*sizeof(double));
    double *Smatrix     = (double *)GpuMallocManaged((size_t)numst * (size_t)numst*sizeof(double));
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

    //    double *vh_x = new double[FP0_BASIS];
    //    double *vh_y = new double[FP0_BASIS];
    //    double *vh_z = new double[FP0_BASIS];
    double *xpsi =(double *)GpuMallocManaged((size_t)numst * (size_t)P0_BASIS*sizeof(double));

    double dipole_ele[3];


    if(0)
    {
        XyzMatrix(Kptr[0], (OrbitalType *)Hmatrix_0, 1, 0, 0);
        if(pct.gridpe == 0)
            for(int i = 0; i < 5; i++) 
            { printf("xyz\n");
                for(int j = 0; j < 5; j++) printf(" %10.4e", 0.001* Hmatrix_0[i*numst + j]);
            }
    }

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
    double time_step = ct.tddft_time_step;

    if(ct.restart_tddft)
    {

        ReadData (ct.infile, vh, rho, vxc, Kptr);
        ReadData_rmgtddft(ct.infile_tddft, vh, vxc, vh_corr, Pn0, Hmatrix, Smatrix,Smatrix, Hmatrix_m1, Hmatrix_0, &pre_steps);
        dcopy(&n2, Hmatrix, &ione, Hmatrix_old, &ione);

    }
    else
    {
        for (int idx = 0; idx < FP0_BASIS; idx++) vtot[idx] = 0.0;
        if(ct.tddft_mode == EFIELD)
        {
            init_efield(vtot);
            GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);
        }
        else if(ct.tddft_mode == POINT_CHARGE)
        {
            init_point_charge_pot(vtot_psi, 1);
        }
        else
        {
            throw RmgFatalException() << " TDDFT mode not defined: "<< ct.tddft_mode<< __FILE__ << " at line " << __LINE__ << "\n";
        }

        HmatrixUpdate(Kptr[0], vtot_psi, (OrbitalType *)matrix_glob);
        Sp->CopySquareMatrixToDistArray(matrix_glob, Akick, numst, desca);

        /* save old vhxc + vnuc */
        for (int idx = 0; idx < FP0_BASIS; idx++) {
            vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];
        }

        // Transfer vtot from the fine grid to the wavefunction grid
        GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);

        /*Generate the Dnm_I */
        get_ddd (vtot, vxc, true);

        HSmatrix (Kptr[0], vtot_psi, vxc_psi, (OrbitalType *)matrix_glob, (OrbitalType *)Smatrix);
        Sp->CopySquareMatrixToDistArray(matrix_glob, Hmatrix, numst, desca);

        my_sync_device();
        dcopy_driver(n2, Hmatrix, ione, Hmatrix_old, ione);

        if(pct.gridpe == 0 && ct.verbose)
        { 
            printf("\nHMa\n");
            for(i = 0; i < 10; i++) 
            {

                printf("\n");
                for(int j = 0; j < 10; j++) printf(" %10.3e",  Hmatrix[i*numst + j]);
            }

            printf("\nSMa\n");
            for(i = 0; i < 10; i++) 
            {

                printf("\n");
                for(int j = 0; j < 10; j++) printf(" %10.3e",  Smatrix[i*numst + j]);
            }
        }


        pre_steps = 0;

        double alpha = 1.0/time_step;
        daxpy_driver ( n2 ,  alpha, Akick, ione , Hmatrix ,  ione) ;

        for(i = 0; i < 2* n2; i++) Pn0[i] = 0.0;

        for(i = 0; i < numst * numst; i++) matrix_glob[i] = 0.0;
        for(i = 0; i < ct.nel/2; i++)  matrix_glob[i * numst + i] = 2.0;
        Sp->CopySquareMatrixToDistArray(matrix_glob, Pn0, numst, desca);


        get_dipole(rho, dipole_ele);

        rmg_printf("\n  x dipolll  %f ", dipole_ele[0]);
        rmg_printf("\n  y dipolll  %f ", dipole_ele[1]);
        rmg_printf("\n  z dipolll  %f ", dipole_ele[2]);

        //   if(pct.gridpe == 0)
        //   for(int i = 0; i < 5; i++) 
        //   { printf("Akick\n");
        //       for(int j = 0; j < 5; j++) printf(" %10.4e", Akick[i*numst + j]);
        //   }


        for(i = 0; i < n2; i++) Hmatrix_m1[i] = 0.0;
        for(i = 0; i < n2; i++) Hmatrix_0[i] = 0.0;
        dcopy_driver(n2, Hmatrix, ione, Hmatrix_m1, ione);
        dcopy_driver(n2, Hmatrix, ione, Hmatrix_0 , ione);
        my_sync_device();
    }

    //  initialize   data for rt-td-dft
    //int nblock = 10 ;   //  size of tthe block for printing (debug!)

    /*
       if(pct.gridpe == 0) { printf("**** Smat  : \n");  print_matrix_d(Smatrix,   &nblock, &numst)   ; }
       if(pct.gridpe == 0) { printf("**** Hmat  : \n");  print_matrix_d(Hmatrix,   &nblock, &numst)   ; }
       if(pct.gridpe == 0) { printf("**** Hmat0 : \n");  print_matrix_d(Hmatrix_0, &nblock, &numst)   ; }
       if(pct.gridpe == 0) { printf("**** Hmat1 : \n");  print_matrix_d(Hmatrix_1, &nblock, &numst)   ; }
     */

    //  run rt-td-dft
    for(tddft_steps = 0; tddft_steps < ct.tddft_steps; tddft_steps++)
    {
        //if(pct.gridpe == 0) printf("=========================================================================\n   step:  %d\n", tddft_steps);

        tot_steps = pre_steps + tddft_steps;

        /* 
        //Wenchang: 
        dscal(&n2, &time_step, Hmatrix, &ione);   
        eldyn_(&numst, Smatrix, Hmatrix, Pn0, Pn1, &Ieldyn, &iprint);
         */
        //  guess H1 from  H(0) and H(-1):
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

            // if(pct.gridpe == 0) { printf("**** Pn1 : \n");   print_matrix_z(Pn1,  &nblock, &numst)  ; }

            //            for(i = 0; i < 10; i++) 
            //            { printf("Pn\n");
            //           for(int j = 0; j < 10; j++) printf(" %8.1e", i, Pn1[i*numst + j]);
            //          }


            /////// <----- update Hamiltonian from  Pn1
            RT2a = new RmgTimer("2-TDDFT: Rho");
            if( scalapack_groups != pct.grid_npes)
            {
                if(Sp->Participates()) 
                {
                    Sp->CopyDistArrayToSquareMatrix(matrix_glob, Pn1, numst, desca);
                    if( scalapack_groups != pct.grid_npes)
                        Sp->ScalapackBlockAllreduce(matrix_glob, (size_t)numst * (size_t)numst);
                }

                Sp->BcastRoot(matrix_glob, numst * numst, MPI_DOUBLE);
            }
            else
            {
                dcopy_driver(n2, Pn1, ione, matrix_glob, ione);
            }


            my_sync_device();

            GetNewRho_rmgtddft((double *)Kptr[0]->orbital_storage, xpsi, rho, matrix_glob, numst);
            delete(RT2a);

            dcopy(&FP0_BASIS, vh_corr, &ione, vh_corr_old, &ione);
            dcopy(&FP0_BASIS, vh, &ione, vh_old, &ione);
            dcopy(&FP0_BASIS, vxc, &ione, vxc_old, &ione);

            //get_vxc(rho, rho_oppo, rhocore, vxc);
            double vtxc, etxc;
            RmgTimer *RT1 = new RmgTimer("2-TDDFT: exchange/correlation");
            Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
            F->v_xc(rho, rhocore, etxc, vtxc, vxc, ct.nspin);
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
            HmatrixUpdate(Kptr[0], vtot_psi, (OrbitalType *)matrix_glob);                                     
            if( scalapack_groups != pct.grid_npes)
            {
                Sp->CopySquareMatrixToDistArray(matrix_glob, Hmatrix, numst, desca);
            }
            else
            {
                my_sync_device();
                dcopy_driver(n2, matrix_glob, ione, Hmatrix, ione);
            }
            delete(RT2a);

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
            RT2a = new RmgTimer("2-TDDFT: Write");
            my_sync_device();
            WriteData_rmgtddft(ct.outfile_tddft, vh, vxc, vh_corr, Pn0, Hmatrix, Smatrix, Smatrix, 
                    Hmatrix_m1, Hmatrix_0, tot_steps+1);
            delete RT2a;
            if(pct.gridpe == 0)fflush(dfi);
        }


    }

    if(pct.gridpe == 0) fclose(dfi);


    RmgTimer *RT2a = new RmgTimer("2-TDDFT: Write");
    my_sync_device();
    WriteData_rmgtddft(ct.outfile_tddft, vh, vxc, vh_corr, Pn0, Hmatrix, Smatrix, Smatrix, 
            Hmatrix_m1, Hmatrix_0, tot_steps+1);
    delete RT2a;
    delete RT0;
}



