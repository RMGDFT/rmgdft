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
#include "RmgSumAll.h"


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

#if USE_NCCL
    int nlocal_ranks;
    MPI_Comm_size(pct.local_comm, &nlocal_ranks);
    if(pct.local_rank == 0)
    {
        ncclGetUniqueId(&ct.nccl_nd_id);
    }
    MPI_Bcast(&ct.nccl_nd_id, sizeof(ct.nccl_nd_id), MPI_BYTE, 0, pct.local_comm);
    ncclCommInitRank(&ct.nccl_local_comm, nlocal_ranks, ct.nccl_nd_id, pct.local_rank);
#endif  

    double *vtot, *vtot_psi;

    int dimx = Rmg_G->get_PX0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimy = Rmg_G->get_PY0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimz = Rmg_G->get_PZ0_GRID(Rmg_G->get_default_FG_RATIO());
    int FP0_BASIS = dimx * dimy * dimz;

    FILE *dfi = NULL, *efi = NULL;
    double vel = get_vel_f();
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
            if(ct.scalapack_block_factor >= ct.num_states - ct.tddft_start_state)
                scalapack_groups = pct.grid_npes;
            break;
        case SUBDIAG_CUSOLVER:
            scalapack_groups = pct.grid_npes;
            break;
        default:
            rmg_error_handler(__FILE__, __LINE__, "Invalid subdiag_driver type in TDDFT");

    } // end switch

#if CUDA_ENABLED || HIP_ENABLED
    if(ct.tddft_gpu)
    {
        scalapack_groups = pct.grid_npes;
    }
    // all tddft propergating will use 1 gpu only, not sure the speed comparison with scalapack for a large system 
#endif
    int last = 1;
    numst = ct.num_states - ct.tddft_start_state; 
    Scalapack *Sp = new Scalapack(scalapack_groups, pct.thisimg, ct.images_per_node, numst,
            ct.scalapack_block_factor, last, pct.grid_comm);
    int Mdim = Sp->GetDistMdim();
    int Ndim = Sp->GetDistNdim();
    n2 = Sp->GetDistMdim() * Sp->GetDistNdim();

    if(!Sp->Participates()) n2 = 1;
    int *desca = Sp->GetDistDesca();

    n22 = 2* n2;

    double *matrix_glob = (double *)RmgMallocHost((size_t)numst * (size_t)numst*sizeof(double));
    double *Hmatrix     = (double *)RmgMallocHost((size_t)n2*sizeof(double));
    double *Hmatrix_old = (double *)RmgMallocHost((size_t)n2*sizeof(double));
    double *Akick       = (double *)RmgMallocHost((size_t)n2*sizeof(double));
    double *Pn0         = (double *)RmgMallocHost((size_t)n2*sizeof(double)*2);
    double *Pn1         = (double *)RmgMallocHost((size_t)n2*sizeof(double)*2);
    double *vh_old      = new double[FP0_BASIS];
    double *vxc_old     = new double[FP0_BASIS];
    double *vh_dipole_old = new double[FP0_BASIS];
    double *vh_dipole     = new double[FP0_BASIS];
    // Jacek: 
    //double *dHmatrix    = new double[n2];   // storage for  H1 -H1_old 
    double *Hmatrix_m1  = (double *)RmgMallocHost((size_t)n2*sizeof(double));
    double *Hmatrix_0   = (double *)RmgMallocHost((size_t)n2*sizeof(double));
    double *Hmatrix_1   = (double *)RmgMallocHost((size_t)n2*sizeof(double));
    double *Hmatrix_dt  = (double *)RmgMallocHost((size_t)n2*sizeof(double));
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

    double dipole_ele[3];

    get_dipole(rho, dipole_ele);
    DipoleCorrection(dipole_ele,  vh_dipole);

    double *Hcore_tddft = (double *)RmgMallocHost((size_t)numst * (size_t)numst*sizeof(double));
    double ES = 0.0, E_downfold = 0.0, EkinPseudo = 0.0, totalE=0.0;
    double ES_0 = 0.0, EkinPseudo_0 = 0.0, totalE_0=0.0;
    double ES_00 = 0.0, EkinPseudo_00 = 0.0, totalE_00=0.0;
    double vtxc, etxc, etxc_0=0.0, etxc_00=0.0;
    std::vector<double> Eterms(6, 1.0);
    double efactor = ct.energy_output_conversion[ct.energy_output_units];
    const char *eunits = ct.energy_output_string[ct.energy_output_units].c_str();

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

    if(pct.gridpe == 0)
    {
        filename = std::string(ct.basename) + "_dipole.dat";

        dfi = fopen(filename.c_str(), "w");

        fprintf(dfi, "\n  &&electric field:  %e  %e  %e ",ct.efield_tddft[0], ct.efield_tddft[1], ct.efield_tddft[2]);
        filename = std::string(ct.basename) + "_totalE";
        efi = fopen(filename.c_str(), "w");

    }




    /* allocate memory for eigenvalue send array and receive array */

    vtot = new double[FP0_BASIS];
    vtot_psi = new double[P0_BASIS];
    fgobj<double> rho_ground;
    double time_step = ct.tddft_time_step;

    int pbasis = Kptr[0]->pbasis;
    size_t psi_alloc = (size_t)ct.num_states * (size_t)pbasis * sizeof(OrbitalType);
    ReadData (ct.infile, vh, rho_ground.data(), vxc, Kptr);
#if CUDA_ENABLED || HIP_ENABLED
    // Wavefunctions are unchanged through TDDFT loop so leave a copy on the GPUs for efficiency.
    // We also need an array of the same size for workspace in HmatrixUpdate and GetNewRho
    gpuMalloc((void **)&Kptr[0]->psi_dev, psi_alloc);
    gpuMalloc((void **)&Kptr[0]->work_dev, psi_alloc);
    gpuMemcpy(Kptr[0]->psi_dev, Kptr[0]->orbital_storage, psi_alloc, gpuMemcpyHostToDevice);
#else
    Kptr[0]->work_cpu = new OrbitalType[psi_alloc];
#endif

    if(ct.restart_tddft)
    {

        ReadData_rmgtddft(ct.infile_tddft, vh, vxc, vh_dipole, Pn0, Hmatrix, Hmatrix_m1, Hmatrix_0, 
                &pre_steps, n2, Eterms, Hcore_tddft, numst);
        E_downfold   = Eterms[0];
        EkinPseudo_0 = Eterms[1];
        ES_0         = Eterms[2];
        etxc_0       = Eterms[3];
        ct.II        = Eterms[4];
        totalE_0     = Eterms[5];
        dcopy(&n2, Hmatrix, &ione, Hmatrix_old, &ione);

    }
    else
    {
        {
            double *v_psi, *vxc_psi;
            int pbasis = Kptr[0]->pbasis;
            v_psi = new double[pbasis];
            vxc_psi = new double[pbasis]();
            int nstates = Kptr[0]->nstates;
            OrbitalType *Hcore = (OrbitalType *)RmgMallocHost(ct.num_kpts_pe * nstates * nstates * sizeof(OrbitalType));

            GetVtotPsi (v_psi, vnuc, Rmg_G->default_FG_RATIO);
            for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
                Kptr[kpt]->ComputeHcore(v_psi, vxc_psi, &Hcore[kpt*nstates * nstates], NULL, NULL);

                for (int st = 0; st < ct.tddft_start_state; st++)
                {
                    double occ = Kptr[kpt]->Kstates[st].occupation[0] * Kptr[kpt]->kp.kweight;
                    E_downfold += std::real(Hcore[kpt * nstates * nstates + st * nstates + st]) * occ;
                }

            }

            for(int st1 = 0; st1 < numst; st1++) {
                for(int st2 = 0; st2 < numst; st2++) {
                    Hcore_tddft[st1 * numst + st2] = std::real(Hcore[(st1+ct.tddft_start_state) * nstates + st2 + ct.tddft_start_state]);
                }
            }

            delete [] v_psi;
            delete [] vxc_psi;
            RmgFreeHost(Hcore);

            //get_vxc(rho, rho_oppo, rhocore, vxc);
            RmgTimer *RT1 = new RmgTimer("2-TDDFT: exchange/correlation");
            Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
            F->v_xc(rho_ground.data(), rhocore, etxc, vtxc, vxc, ct.nspin);
            etxc_00 = etxc;
            delete F;
            delete RT1;

            RT1 = new RmgTimer("2-TDDFT: Vh");
            VhDriver(rho_ground.data(), rhoc, vh, ct.vh_ext, 1.0-12);
            delete RT1;

            ES_00 = 0.0;
            for (int idx = 0; idx < FP0_BASIS; idx++) 
            {
                ES_00 += (rho_ground[idx] - rhoc[idx]) * vh[idx];
            }

            ES_00 = 0.5 * vel * RmgSumAll(ES_00, pct.grid_comm);
            int ntot2 = numst *numst;
            int ione = 1;
            for(i = 0; i < numst * numst; i++) matrix_glob[i] = 0.0;
            for(int i = 0; i < numst; i++) matrix_glob[i * numst + i] = Kptr[0]->Kstates[i + ct.tddft_start_state].occupation[0];
            EkinPseudo_00 = ddot(&ntot2, matrix_glob, &ione, Hcore_tddft, &ione);
            totalE_00 = E_downfold + EkinPseudo_00 + ES_00 + etxc_00 + ct.II;

        }
        for (int idx = 0; idx < FP0_BASIS; idx++) vtot[idx] = 0.0;
        if(ct.tddft_mode == EFIELD)
        {
            init_efield(vtot, ct.efield_tddft);
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


        HmatrixUpdate(Kptr[0], vtot_psi, (OrbitalType *)matrix_glob, ct.tddft_start_state);
        Sp->CopySquareMatrixToDistArray(matrix_glob, Akick, numst, desca);


        for(int i = 0; i < numst * numst; i++) matrix_glob[i] = 0.0; 
        for(int i = 0; i < numst; i++) matrix_glob[i * numst + i] = Kptr[0]->Kstates[i + ct.tddft_start_state].eig[0];

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

        }


        pre_steps = 0;

        double alpha = 1.0/time_step;
        daxpy_driver ( n2 ,  alpha, Akick, ione , Hmatrix ,  ione) ;

        for(i = 0; i < 2* n2; i++) Pn0[i] = 0.0;

        for(i = 0; i < numst * numst; i++) matrix_glob[i] = 0.0;
        for(int i = 0; i < numst; i++) matrix_glob[i * numst + i] = Kptr[0]->Kstates[i + ct.tddft_start_state].occupation[0];
        Sp->CopySquareMatrixToDistArray(matrix_glob, Pn0, numst, desca);


        if(pct.gridpe == 0)
        {
            fprintf(dfi, "\n  &&dipole at groud state:  %18.10e  %18.10e  %18.10e ",
                dipole_ele[0], dipole_ele[1], dipole_ele[2]);
        }
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
       if(pct.gridpe == 0) { printf("**** Hmat  : \n");  print_matrix_d(Hmatrix,   &nblock, &numst)   ; }
       if(pct.gridpe == 0) { printf("**** Hmat0 : \n");  print_matrix_d(Hmatrix_0, &nblock, &numst)   ; }
       if(pct.gridpe == 0) { printf("**** Hmat1 : \n");  print_matrix_d(Hmatrix_1, &nblock, &numst)   ; }
     */

    //  run rt-td-dft
    for(tddft_steps = 0; tddft_steps < ct.tddft_steps; tddft_steps++)
    {
        //if(pct.gridpe == 0) printf("=========================================================================\n   step:  %d\n", tddft_steps);

        tot_steps = pre_steps + tddft_steps;

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
            GetNewRho_rmgtddft(Kptr[0], rho, matrix_glob, numst, ct.tddft_start_state, rho_ground.data());
            delete(RT2a);

            dcopy(&FP0_BASIS, vh_dipole, &ione, vh_dipole_old, &ione);
            dcopy(&FP0_BASIS, vh, &ione, vh_old, &ione);
            dcopy(&FP0_BASIS, vxc, &ione, vxc_old, &ione);

            //get_vxc(rho, rho_oppo, rhocore, vxc);
            RmgTimer *RT1 = new RmgTimer("2-TDDFT: exchange/correlation");
            Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
            F->v_xc(rho, rhocore, etxc, vtxc, vxc, ct.nspin);
            delete F;
            delete RT1;

            RT1 = new RmgTimer("2-TDDFT: Vh");
            VhDriver(rho, rhoc, vh, ct.vh_ext, 1.0-12);
            delete RT1;

            get_dipole(rho, dipole_ele);
            DipoleCorrection(dipole_ele,  vh_dipole);

            for (int idx = 0; idx < FP0_BASIS; idx++) {
                vtot[idx] = vxc[idx] + vh[idx] + vh_dipole[idx]
                    -vxc_old[idx] -vh_old[idx] - vh_dipole_old[idx];
            }

            ES = 0.0;
            for (int idx = 0; idx < FP0_BASIS; idx++) 
            {
                ES += (rho[idx] - rhoc[idx]) * vh[idx];
            }

            ES = 0.5 * vel * RmgSumAll(ES, pct.grid_comm);
            int ntot2 = numst *numst;
            int ione = 1;
            EkinPseudo = ddot(&ntot2, matrix_glob, &ione, Hcore_tddft, &ione);
            totalE = E_downfold + EkinPseudo + ES + etxc + ct.II;

            if(tot_steps == 0 )
            {
                totalE_0 = totalE;
                EkinPseudo_0 = EkinPseudo;
                ES_0 = ES;
                etxc_0 = etxc;
            }


            GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);
            RT2a = new RmgTimer("2-TDDFT: Hupdate");
            HmatrixUpdate(Kptr[0], vtot_psi, (OrbitalType *)matrix_glob, ct.tddft_start_state);                                     
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
            rmg_printf("step: %5d  iteration: %d  thrs= %12.5e err=  %12.5e at element: %5d \n", 
                    tddft_steps, iter_scf,    thrs_dHmat,  err,         ij_err);  
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

        if(pct.gridpe == 0)
        {
            if(tot_steps == 0) 
            {
                Eterms[0] = E_downfold  ;
                Eterms[1] = EkinPseudo_0;
                Eterms[2] = ES_0        ;
                Eterms[3] = etxc_0      ;
                Eterms[4] = ct.II       ;
                Eterms[5] = totalE_0    ;
                fprintf(efi, " && totalE_0, EkinPseudo_0, Vh_0, Exc_0  %s", eunits);
                fprintf(efi, "\n&& %18.10f  %18.10f  %18.10f  %18.10f at Ground state", totalE_00, EkinPseudo_00, ES_00, etxc_00);
                fprintf(efi, "\n&& %18.10f  %18.10f  %18.10f  %18.10f at 1st TDDFT step", totalE_0, EkinPseudo_0, ES_0, etxc_0);
                fprintf(efi, "\n&&time  EkinPseudo_diff, Vh_diff, Exc_diff  totalE_diff%s", eunits);
            }
            fprintf(efi, "\n  %f  %16.8e %16.8e,%16.8e,%16.8e   ",
                    tot_steps*time_step, (EkinPseudo-EkinPseudo_0) * efactor, (ES-ES_0) * efactor, (etxc-etxc_0) * efactor, (totalE-totalE_0) * efactor);
        }
        if(pct.gridpe == 0)fprintf(dfi, "\n  %f  %18.10e  %18.10e  %18.10e ",
                tot_steps*time_step, dipole_ele[0], dipole_ele[1], dipole_ele[2]);

        if((tddft_steps +1) % ct.checkpoint == 0)
        {   
            RT2a = new RmgTimer("2-TDDFT: Write");
            my_sync_device();
            WriteData_rmgtddft(ct.outfile_tddft, vh, vxc, vh_dipole, Pn0, Hmatrix, 
                    Hmatrix_m1, Hmatrix_0, tot_steps+1, n2, Eterms, Hcore_tddft, numst);
            delete RT2a;
            if(pct.gridpe == 0)fflush(dfi);
            if(pct.gridpe == 0)fflush(efi);
        }


    }

#if CUDA_ENABLED || HIP_ENABLED
    gpuFree(Kptr[0]->work_dev);
    gpuFree(Kptr[0]->psi_dev);
#endif
    if(pct.gridpe == 0) 
    {
        fclose(dfi);
        fclose(efi);
    }

    RmgTimer *RT2a = new RmgTimer("2-TDDFT: Write");
    my_sync_device();
    WriteData_rmgtddft(ct.outfile_tddft, vh, vxc, vh_dipole, Pn0, Hmatrix, 
            Hmatrix_m1, Hmatrix_0, tot_steps+1, n2, Eterms, Hcore_tddft, numst);
    delete RT2a;
    delete RT0;
}


