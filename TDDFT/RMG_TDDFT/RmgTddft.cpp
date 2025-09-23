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
void eldyn_ort(int *desca, int Mdim, int Ndim, std::complex<double> *F,std::complex<double> *Po0,std::complex<double> *Po1,int *p_Ieldyn,  double *thrs,int*maxiter,  double *errmax,int
        *niter , int *p_iprint, MPI_Comm comm) ;
void eldyn_nonort(int *p_N, double *S, double *F,double *Po0,double *Pn1,int *p_Ieldyn,  double *thrs,int*maxiter,  double *errmax,int *niter , int *p_iprint) ;
void  tstconv(double *C,int *p_M, double *p_thrs,int *p_ierr, double *p_err, bool *p_tconv, MPI_Comm comm);

void CalculatePrho(double *, double *);

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


template void RmgTddft<double> ( spinobj<double> &vxc,
                   fgobj<double> &vh,
                   fgobj<double> &vnuc,
                   spinobj<double> &rho,
                   fgobj<double> &rhocore,
                   fgobj<double> &rhoc,Kpoint<double> **Kptr);

template void RmgTddft<std::complex<double> > ( spinobj<double> &vxc,
                   fgobj<double> &vh,
                   fgobj<double> &vnuc,
                   spinobj<double> &rho,
                   fgobj<double> &rhocore,
                   fgobj<double> &rhoc, Kpoint<std::complex<double>> **);
template <typename OrbitalType> void RmgTddft ( spinobj<double> &vxc,
                   fgobj<double> &vh,
                   fgobj<double> &vnuc,
                   spinobj<double> &rho,
                   fgobj<double> &rhocore,
                   fgobj<double> &rhoc,Kpoint<OrbitalType> **Kptr)
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

    double dipole_tot[3];

    int dimx = Rmg_G->get_PX0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimy = Rmg_G->get_PY0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimz = Rmg_G->get_PZ0_GRID(Rmg_G->get_default_FG_RATIO());
    int FP0_BASIS = dimx * dimy * dimz;

    FILE *dfi = NULL, *efi = NULL, *current_fi = NULL, *dbp_fi = NULL;
    double vel = get_vel_f();
    std::string filename;
    int n2,n22, numst, P0_BASIS,i, ione =1;
    int tot_steps = 0, pre_steps = 0, tddft_steps;
    int Ieldyn = 1;    // BCH  
                       //int Ieldyn = 2;    // Diagev
    int iprint = ct.verbose;



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

    if(ct.tddft_gpu)
    {
        scalapack_groups = pct.grid_npes;
    }
    // all tddft propergating will use 1 gpu only, not sure the speed comparison with scalapack for a large system 
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
    int n2_C = n2 * sizeof(OrbitalType)/sizeof(double);
    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        Kptr[kpt]->Hmatrix_cpu     = (OrbitalType *)RmgMallocHost((size_t)n2*sizeof(OrbitalType));
        Kptr[kpt]->Pn0_cpu         = (OrbitalType *)RmgMallocHost((size_t)n2*sizeof(double)*2);
        Kptr[kpt]->Pn1_cpu         = (OrbitalType *)RmgMallocHost((size_t)n2*sizeof(double)*2);
        Kptr[kpt]->Hmatrix_m1_cpu  = (OrbitalType *)RmgMallocHost((size_t)n2*sizeof(OrbitalType));
        Kptr[kpt]->Hmatrix_1_cpu  = (OrbitalType *)RmgMallocHost((size_t)n2*sizeof(OrbitalType));
        Kptr[kpt]->Hmatrix_0_cpu   = (OrbitalType *)RmgMallocHost((size_t)n2*sizeof(OrbitalType));
        if(ct.tddft_mode == VECTOR_POT)
        {
            Kptr[kpt]->Pxmatrix_cpu   = (OrbitalType *)RmgMallocHost((size_t)n2*sizeof(OrbitalType));
            Kptr[kpt]->Pymatrix_cpu   = (OrbitalType *)RmgMallocHost((size_t)n2*sizeof(OrbitalType));
            Kptr[kpt]->Pzmatrix_cpu   = (OrbitalType *)RmgMallocHost((size_t)n2*sizeof(OrbitalType));

        }
        else
        {
            Kptr[kpt]->Akick_cpu   = (OrbitalType *)RmgMallocHost((size_t)n2*sizeof(OrbitalType));
        }
    }

    OrbitalType *matrix_glob = (OrbitalType *)RmgMallocHost((size_t)numst * (size_t)numst*sizeof(OrbitalType));
    fgobj<double> vh_old, vh_dipole, vh_dipole_old;
    spinobj<double> vxc_old, rho_k, rho_ksum;
    // Jacek: 
    //double *dHmatrix    = new double[n2];   // storage for  H1 -H1_old 
    OrbitalType *Pn1        ;
    OrbitalType *Hmatrix_1  ;

    OrbitalType *Hmatrix     = NULL;
    OrbitalType *Pn0         = NULL;
    OrbitalType *Hmatrix_m1  = NULL;
    OrbitalType *Hmatrix_0   = NULL;
    size_t matrix_size = n2*sizeof(OrbitalType);
    if(ct.tddft_gpu)
    {
        gpuMalloc((void **)&Hmatrix, n2*sizeof(OrbitalType));
        gpuMalloc((void **)&Hmatrix_m1, n2*sizeof(OrbitalType));
        gpuMalloc((void **)&Hmatrix_0, n2*sizeof(OrbitalType));
        gpuMalloc((void **)&Pn0, 2*n2*sizeof(double));

        gpuMalloc((void **)&Hmatrix_1, n2*sizeof(OrbitalType));
        gpuMalloc((void **)&Pn1, 2*n2*sizeof(double));
    }
    else
    {
        Pn1  = (OrbitalType *)RmgMallocHost((size_t)n2*sizeof(double)*2);
    }


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



    get_dipole(rho.data(), rhoc.data(), dipole_tot);
    if(ct.dipole_corr[0]+ct.dipole_corr[1]+ct.dipole_corr[2] >0)
    {
        DipoleCorrection(dipole_tot,  vh_dipole.data());
    }

    double *Hcore_tddft = new double[(size_t)numst * (size_t)numst];
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

    if(!ct.norm_conserving_pp)
    {
        rmg_printf(" \n  TDDFT support NCPP only \n");

        fflush(NULL);
        throw RmgFatalException() << " TDDFT support NCPP only in "<< __FILE__ << " at line " << __LINE__ << "\n";
        exit(0);
    }

    if(pct.kstart == 0 && pct.gridpe == 0)
    {

        
        if(ct.tddft_mode == VECTOR_POT)
        {
            filename = std::string(ct.basename)+"_spin" +std::to_string(pct.spinpe)+ "_current.dat";
            
            current_fi = fopen(filename.c_str(), "w");
            fprintf(current_fi, "\n  &&electric field in cartesian unit:  %e  %e  %e ",ct.efield_tddft_crds[0], ct.efield_tddft_crds[1], ct.efield_tddft_crds[2]);

            if(ct.BerryPhase)
            {
                filename = std::string(ct.basename) +"_spin" +std::to_string(pct.spinpe)+ "_bp_dipole.dat";
                dbp_fi = fopen(filename.c_str(), "w");
                fprintf(dbp_fi, "\n  &&electric field in cartesian unit:  %e  %e  %e ",ct.efield_tddft_crds[0], ct.efield_tddft_crds[1], ct.efield_tddft_crds[2]);
            }
        }
        else
        {
            filename = std::string(ct.basename) +"_spin" +std::to_string(pct.spinpe)+ "_dipole.dat";

            dfi = fopen(filename.c_str(), "w");

            fprintf(dfi, "\n  &&electric field in cartesian unit:  %e  %e  %e ",ct.efield_tddft_crds[0], ct.efield_tddft_crds[1], ct.efield_tddft_crds[2]);
            filename = std::string(ct.basename) + "_totalE";
            efi = fopen(filename.c_str(), "w");
        }
    }



    /* allocate memory for eigenvalue send array and receive array */

    wfobj<double> vtot_psi;
    fgobj<double> vtot;
    spinobj<double> rho_ground;
    double time_step = ct.tddft_time_step;

    int pbasis = Kptr[0]->pbasis;
    size_t psi_alloc = (size_t)ct.num_states * (size_t)pbasis * sizeof(OrbitalType);
    ReadData (ct.infile, vh.data(), rho_ground.data(), vxc.data(), Kptr);
    rho_ground.get_oppo();

    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
#if CUDA_ENABLED || HIP_ENABLED
        // Wavefunctions are unchanged through TDDFT loop so leave a copy on the GPUs for efficiency.
        // We also need an array of the same size for workspace in HmatrixUpdate and GetNewRho
        gpuMalloc((void **)&Kptr[kpt]->psi_dev, psi_alloc);
        gpuMalloc((void **)&Kptr[kpt]->work_dev, psi_alloc);
        RmgMemcpy(Kptr[kpt]->psi_dev, Kptr[kpt]->orbital_storage, psi_alloc);
#else
        Kptr[kpt]->work_cpu = new OrbitalType[(size_t)ct.num_states * (size_t)pbasis];
#endif
    }
    if(ct.restart_tddft)
    {

        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {
            int kpt_glob = kpt + pct.kstart;

            char newname[MAX_PATH + 20];
            sprintf (newname, "%s_spin%d_kpt%d_gridpe%d", ct.infile_tddft, pct.spinpe, kpt_glob, pct.gridpe);
            ReadData_rmgtddft(newname, vh.data(), vxc.data(), vh_dipole.data(), (double *)Kptr[kpt]->Pn0_cpu, (double *)Kptr[kpt]->Hmatrix_cpu, 
                    (double *)Kptr[kpt]->Hmatrix_m1_cpu, (double *)Kptr[kpt]->Hmatrix_0_cpu, 
                    &pre_steps, n2, n2_C, Eterms, Hcore_tddft, numst);
        }
        E_downfold   = Eterms[0];
        EkinPseudo_0 = Eterms[1];
        ES_0         = Eterms[2];
        etxc_0       = Eterms[3];
        ct.II        = Eterms[4];
        totalE_0     = Eterms[5];
    }
    else
    {
        {
            int pbasis = Kptr[0]->pbasis;
            wfobj<double> v_psi, vxc_psi;
            int nstates = Kptr[0]->nstates;
            OrbitalType *Hcore = (OrbitalType *)RmgMallocHost(ct.num_kpts_pe * nstates * nstates * sizeof(OrbitalType));

            GetVtotPsi (v_psi.data(), vnuc.data(), Rmg_G->default_FG_RATIO);
            for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
                Kptr[kpt]->ComputeHcore(v_psi.data(), vxc_psi.data(), &Hcore[kpt*nstates * nstates], NULL, NULL);

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

            RmgFreeHost(Hcore);

            //get_vxc(rho, rho_oppo, rhocore, vxc);
            RmgTimer *RT1 = new RmgTimer("2-TDDFT: exchange/correlation");
            compute_vxc(rho_ground.data(), rhocore.data(), etxc, vtxc, vxc.data(), ct.nspin);
            etxc_00 = etxc;
            delete RT1;

            RT1 = new RmgTimer("2-TDDFT: Vh");
            VhDriver(rho_ground.data(), rhoc.data(), vh.data(), ct.vh_ext, 1.0-12);
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
            EkinPseudo_00 = ddot(&ntot2, (double *)matrix_glob, &ione, Hcore_tddft, &ione);
            totalE_00 = E_downfold + EkinPseudo_00 + ES_00 + etxc_00 + ct.II;

        }

        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {

            for(int i = 0; i < numst * numst; i++) matrix_glob[i] = 0.0; 
            for(int i = 0; i < numst; i++) matrix_glob[i * numst + i] = Kptr[kpt]->Kstates[i + ct.tddft_start_state].eig[0];

            Sp->CopySquareMatrixToDistArray(matrix_glob, Kptr[kpt]->Hmatrix_cpu, numst, desca);
            memcpy(Kptr[kpt]->Hmatrix_1_cpu, Kptr[kpt]->Hmatrix_cpu, matrix_size);
            memcpy(Kptr[kpt]->Hmatrix_0_cpu, Kptr[kpt]->Hmatrix_cpu, matrix_size);
            memcpy(Kptr[kpt]->Hmatrix_m1_cpu, Kptr[kpt]->Hmatrix_0_cpu, matrix_size);


            if(n2 == n2_C)
            {
                for(i = 0; i < 2* n2; i++) Kptr[kpt]->Pn0_cpu[i] = 0.0;
            }
            else
            {
                for(i = 0; i < n2; i++) Kptr[kpt]->Pn0_cpu[i] = 0.0;
            }

            for(i = 0; i < numst * numst; i++) matrix_glob[i] = 0.0;
            for(int i = 0; i < numst; i++) matrix_glob[i * numst + i] = Kptr[kpt]->Kstates[i + ct.tddft_start_state].occupation[0];
            Sp->CopySquareMatrixToDistArray(matrix_glob, Kptr[kpt]->Pn0_cpu, numst, desca);
        }

        rmg_printf("\n  x dipolll  %f ", dipole_tot[0]);
        rmg_printf("\n  y dipolll  %f ", dipole_tot[1]);
        rmg_printf("\n  z dipolll  %f \n", dipole_tot[2]);
        fflush(NULL);

        //   if(pct.gridpe == 0)
        //   for(int i = 0; i < 5; i++) 
        //   { printf("Akick\n");
        //       for(int j = 0; j < 5; j++) printf(" %10.4e", Akick[i*numst + j]);
        //   }

        if(ct.tddft_mode == EFIELD || ct.tddft_mode == POINT_CHARGE)
        {
            double alpha = 1.0/time_step;
            for (int idx = 0; idx < FP0_BASIS; idx++) vtot[idx] = 0.0;
            if(ct.tddft_mode == EFIELD)
            {
                init_efield(vtot.data(), ct.efield_tddft_crds);
                GetVtotPsi (vtot_psi.data(), vtot.data(), Rmg_G->default_FG_RATIO);
            }
            else if(ct.tddft_mode == POINT_CHARGE)
            {
                init_point_charge_pot(vtot_psi.data(), 1);
            }
            for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
                HmatrixUpdate(Kptr[kpt], vtot_psi.data(), (OrbitalType *)matrix_glob, ct.tddft_start_state);
                Sp->CopySquareMatrixToDistArray(matrix_glob, Kptr[kpt]->Akick_cpu, numst, desca);
                daxpy ( &n2_C,  &alpha, (double *)Kptr[kpt]->Akick_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_0_cpu , &ione) ;
                memcpy(Kptr[kpt]->Hmatrix_m1_cpu, Kptr[kpt]->Hmatrix_0_cpu, matrix_size);
            }
        }
    }

    //  initialize   data for rt-td-dft
    //int nblock = 10 ;   //  size of tthe block for printing (debug!)

    /*
       if(pct.gridpe == 0) { printf("**** Hmat  : \n");  print_matrix_d(Hmatrix,   &nblock, &numst)   ; }
       if(pct.gridpe == 0) { printf("**** Hmat0 : \n");  print_matrix_d(Hmatrix_0, &nblock, &numst)   ; }
       if(pct.gridpe == 0) { printf("**** Hmat1 : \n");  print_matrix_d(Hmatrix_1, &nblock, &numst)   ; }
     */

    if(ct.tddft_mode == VECTOR_POT)
    {
        // vector potential will be A(t) =  ct.efield_tddft * cos(tddft_frequency * t)
        // VecP matrix is <psi| ct.efied_tddft dot gradient | psi> 
        //
        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
            VecPHmatrix(Kptr[kpt], ct.efield_tddft_crds, desca, ct.tddft_start_state);
            if(pre_steps == 0)
            {
                // at t= 0, cos(omega t) = 1.0
                daxpy ( &n2_C ,  &ct.efield_tddft_crds[0], (double *)Kptr[kpt]->Pxmatrix_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_m1_cpu,  &ione) ;
                daxpy ( &n2_C ,  &ct.efield_tddft_crds[1], (double *)Kptr[kpt]->Pymatrix_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_m1_cpu,  &ione) ;
                daxpy ( &n2_C ,  &ct.efield_tddft_crds[2], (double *)Kptr[kpt]->Pzmatrix_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_m1_cpu,  &ione) ;
                memcpy(Kptr[kpt]->Hmatrix_0_cpu, Kptr[kpt]->Hmatrix_m1_cpu, matrix_size);
            }

            CurrentNlpp(Kptr[kpt], desca, ct.tddft_start_state);
        }

    }

    double current[3], current0[3];
    current[0] = 0.0;
    current[1] = 0.0;
    current[2] = 0.0;
    current0[0] = 0.0;
    current0[1] = 0.0;
    current0[2] = 0.0;
    double tot_bp_pol = 0.0;
    if(ct.tddft_mode == VECTOR_POT )
    {
        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
            std::complex<double> tem_x = zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->Pxmatrix_cpu, &ione);
            current0[0] += std::real(tem_x) * Kptr[kpt]->kp.kweight;
            std::complex<double> tem_y = zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->Pymatrix_cpu, &ione);
            current0[1] += std::real(tem_y) * Kptr[kpt]->kp.kweight;
            std::complex<double> tem_z = zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->Pzmatrix_cpu, &ione);
            current0[2] += std::real(tem_z) * Kptr[kpt]->kp.kweight;
        }
        if(ct.BerryPhase)
        {
            // Rmg_BP->CalcBP_Skk1(Kptr, ct.tddft_start_state, matrix_glob, *Sp);
            // Rmg_BP->CalcBP_tddft(Kptr, tot_bp_pol, matrix_glob, *Sp);
            Rmg_BP->tddft_Xml(Kptr, ct.tddft_start_state, matrix_glob, *Sp);
            tot_bp_pol = 0.0;
            for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
                std::complex<double> tem_x = zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->BP_Xml, &ione);
                tot_bp_pol += std::real(tem_x) * Kptr[kpt]->kp.kweight;
            }
            MPI_Allreduce(MPI_IN_PLACE, &tot_bp_pol, 1, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
            Sp->ScalapackBlockAllreduce(&tot_bp_pol, 1);
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, current0, 3, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
    Sp->ScalapackBlockAllreduce(current0, 3);
    Rmg_Symm->symm_vec(current0);

    if(pct.kstart == 0 && pct.gridpe == 0)
    {
        if(ct.tddft_mode == VECTOR_POT)
        {
            fprintf(current_fi, "\n  &&current at groud state:  %18.10e  %18.10e  %18.10e nonzero due to kpoint sampling",
                    current0[0], current0[1], current0[2]);
        }
        else
        {
            fprintf(dfi, "\n  &&dipole at groud state:  %18.10e  %18.10e  %18.10e ",
                    dipole_tot[0], dipole_tot[1], dipole_tot[2]);
        }
        if(ct.BerryPhase)
        {
            fprintf(dbp_fi, "\n  &&dipole at groud state BeryPhase (C/m^2):  %18.10e  %18.10e  %18.10e ",
                    tot_bp_pol, 0.0,0.0);
        }
    }
    //  run rt-td-dft
    RmgTimer *RT2a ;    // timer type  declaration
    for(tddft_steps = 0; tddft_steps < ct.tddft_steps; tddft_steps++)
    {
        //if(pct.gridpe == 0) printf("=========================================================================\n   step:  %d\n", tddft_steps);

        tot_steps = pre_steps + tddft_steps;

        //  guess H1 from  H(0) and H(-1):

        current[0] = 0.0;
        current[1] = 0.0;
        current[2] = 0.0;

        RT2a = new RmgTimer("2-TDDFT: extrapolate");
        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
            //if(ct.tddft_mode == VECTOR_POT && tot_steps == 0)
            //{
            //    double coswt = cos(ct.tddft_frequency * tot_steps * time_step);
            //    double coswtx = coswt * ct.efield_tddft_crds[0];
            //    double coswty = coswt * ct.efield_tddft_crds[1];
            //    double coswtz = coswt * ct.efield_tddft_crds[2];
            //    daxpy ( &n2_C ,  &coswtx, (double *)Kptr[kpt]->Pxmatrix_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_0_cpu,  &ione) ;
            //    daxpy ( &n2_C ,  &coswty, (double *)Kptr[kpt]->Pymatrix_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_0_cpu,  &ione) ;
            //    daxpy ( &n2_C ,  &coswtz, (double *)Kptr[kpt]->Pzmatrix_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_0_cpu,  &ione) ;
            //}
            extrapolate_Hmatrix ((double *)Kptr[kpt]->Hmatrix_m1_cpu, (double *)Kptr[kpt]->Hmatrix_0_cpu, (double *)Kptr[kpt]->Hmatrix_1_cpu, n2_C) ;
        }   

        my_sync_device();
        delete RT2a;


        int  Max_iter_scf = 10 ; int  iter_scf =0 ;
        err =1.0e0   ;  thrs_dHmat  = 1e-7  ;

        double  thrs_bch =1.0e-7; 
        int     maxiter_bch  =100;
        double  errmax_bch ;
        int     niter_bch ;


        //-----   SCF loop  starts here: 
        while (err > thrs_dHmat &&  iter_scf <  Max_iter_scf)  {

            for(int idx = 0; idx < FP0_BASIS; idx++) rho_ksum[idx] = 0.0;
            for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
                RT2a = new RmgTimer("2-TDDFT: memcpy");
                if(ct.tddft_gpu)
                {
                    RmgMemcpy(Hmatrix, Kptr[kpt]->Hmatrix_cpu, matrix_size);
                    RmgMemcpy(Hmatrix_m1, Kptr[kpt]->Hmatrix_m1_cpu, matrix_size);
                    RmgMemcpy(Hmatrix_1, Kptr[kpt]->Hmatrix_1_cpu, matrix_size);
                    RmgMemcpy(Hmatrix_0, Kptr[kpt]->Hmatrix_0_cpu, matrix_size);
                    RmgMemcpy(Pn0, Kptr[kpt]->Pn0_cpu, 2*n2*sizeof(double));
                }
                else
                {
                    Hmatrix = Kptr[kpt]->Hmatrix_cpu;
                    Hmatrix_m1 = Kptr[kpt]->Hmatrix_m1_cpu;
                    Hmatrix_1 = Kptr[kpt]->Hmatrix_1_cpu;
                    Hmatrix_0 = Kptr[kpt]->Hmatrix_0_cpu;
                    Pn0 = Kptr[kpt]->Pn0_cpu;
                    Pn1 = Kptr[kpt]->Pn1_cpu;
                }
                my_sync_device();
                delete RT2a;
                RT2a = new RmgTimer("2-TDDFT: ELDYN");
                magnus ((double *)Hmatrix_0,    (double *)Hmatrix_1 , time_step, (double *)Hmatrix_m1 , n2_C) ; 
                /* --- C++  version:  --*/

                eldyn_ort(desca, Mdim, Ndim,  Hmatrix_m1,Pn0,Pn1,&Ieldyn, &thrs_bch,&maxiter_bch,  &errmax_bch,&niter_bch ,  &iprint, Sp->GetComm()) ;
                RmgMemcpy(Kptr[kpt]->Pn1_cpu, Pn1, 2*n2*sizeof(double));

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
                            Sp->ScalapackBlockAllreduce((double *)matrix_glob, (size_t)numst * (size_t)numst *n2_C/n2);
                    }

                    size_t count = numst * numst * sizeof(OrbitalType) /sizeof(double);
                    Sp->BcastRoot((double *)matrix_glob, count, MPI_DOUBLE);
                }
                else
                {
                    RmgMemcpy(matrix_glob, Pn1, matrix_size);
                }


                my_sync_device();
                GetNewRho_rmgtddft(Kptr[kpt], rho_k.data(), matrix_glob, numst, ct.tddft_start_state);

                int kpt_glob = kpt + pct.kstart;
                for(int idx = 0; idx < FP0_BASIS; idx++) rho_ksum[idx] += rho_k[idx] * ct.kp[kpt_glob].kweight;

                delete(RT2a);

            }

            MPI_Allreduce(MPI_IN_PLACE, rho_ksum.data(), FP0_BASIS, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
            for(int idx = 0; idx < FP0_BASIS; idx++) rho[idx] = rho_ksum[idx] + rho_ground[idx];
            rho.get_oppo();


            //write_rho_x(rho, "update rho");
            //write_rho_x(rho_ground.data(), "groumd rho");

            dcopy(&FP0_BASIS, vh_dipole.data(), &ione, vh_dipole_old.data(), &ione);
            dcopy(&FP0_BASIS, vh.data(), &ione, vh_old.data(), &ione);
            dcopy(&FP0_BASIS, vxc.data(), &ione, vxc_old.data(), &ione);

            //get_vxc(rho, rho_oppo, rhocore, vxc);
            RmgTimer *RT1 = new RmgTimer("2-TDDFT: exchange/correlation");
            compute_vxc(rho.data(), rhocore.data(), etxc, vtxc, vxc.data(), ct.nspin);
            delete RT1;

            RT1 = new RmgTimer("2-TDDFT: Vh");
            VhDriver(rho.data(), rhoc.data(), vh.data(), ct.vh_ext, 1.0-12);
            delete RT1;

            get_dipole(rho.data(), rhoc.data(), dipole_tot);
            if(ct.dipole_corr[0]+ct.dipole_corr[1]+ct.dipole_corr[2] >0)
            {
                DipoleCorrection(dipole_tot,  vh_dipole.data());
            }

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
            EkinPseudo = ddot(&ntot2, (double *)matrix_glob, &ione, Hcore_tddft, &ione);
            totalE = E_downfold + EkinPseudo + ES + etxc + ct.II;

            if(tot_steps == 0 )
            {
                totalE_0 = totalE;
                EkinPseudo_0 = EkinPseudo;
                ES_0 = ES;
                etxc_0 = etxc;
            }


            GetVtotPsi (vtot_psi.data(), vtot.data(), Rmg_G->default_FG_RATIO);

            for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
                RT2a = new RmgTimer("2-TDDFT: Hupdate");
                HmatrixUpdate(Kptr[kpt], vtot_psi.data(), matrix_glob, ct.tddft_start_state);                                     
                if( scalapack_groups != pct.grid_npes)
                {
                    Sp->CopySquareMatrixToDistArray(matrix_glob, Kptr[kpt]->Hmatrix_m1_cpu, numst, desca);
                }
                else
                {
                    memcpy(Kptr[kpt]->Hmatrix_m1_cpu, matrix_glob, matrix_size);
                }
                delete(RT2a);

                RT2a = new RmgTimer("2-TDDFT: conv check");
                my_sync_device();
                double one = 1.0, mone = -1.0;
                daxpy( &n2_C ,  &one, (double *)Kptr[kpt]->Hmatrix_m1_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_cpu,  &ione) ;

                //////////  < ---  end of Hamiltonian update

                // check error and update Hmatrix_1:
                my_sync_device();
                daxpy ( &n2_C ,  &mone, (double *)Kptr[kpt]->Hmatrix_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_1_cpu ,  &ione) ;

                //tst_conv_matrix (&err, &ij_err ,  Hmatrix_1,  n2, Sp->GetComm()) ;  //  check error  how close  H and H_old are

                bool tConv;
                tstconv((double *)Kptr[kpt]->Hmatrix_1_cpu, &n2_C, &thrs_dHmat,&ij_err,&err,&tConv, Sp->GetComm());
                memcpy(Kptr[kpt]->Hmatrix_1_cpu, Kptr[kpt]->Hmatrix_cpu, matrix_size);
                delete(RT2a);
            }

            MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, pct.kpsub_comm);


            if(pct.imgpe == 0) { printf("step: %5d  iteration: %d  thrs= %12.5e err=  %12.5e at element: %5d \n", 
                    tddft_steps, iter_scf,    thrs_dHmat,  err,         ij_err); } 
            rmg_printf("step: %5d  iteration: %d  thrs= %12.5e err=  %12.5e at element: %5d \n", 
                    tddft_steps, iter_scf,    thrs_dHmat,  err,         ij_err);  
            //err= -1.0e0 ;  
            iter_scf ++ ;
        } //---- end of  SCF/while loop 


        RT2a = new RmgTimer("2-TDDFT: current and dipole");
        //  extract dipole from rho(Pn1)
        get_dipole(rho.data(), rhoc.data(), dipole_tot);
        /*  done with propagation,  save Pn1 ->  Pn0 */

        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
            memcpy(Kptr[kpt]->Pn0_cpu, Kptr[kpt]->Pn1_cpu, n22 * sizeof(double));

            // save current  H0, H1 for the  next step extrapolatiion
            memcpy(Kptr[kpt]->Hmatrix_m1_cpu, Kptr[kpt]->Hmatrix_0_cpu, matrix_size);
            //dcopy(&n2, Hmatrix  , &ione, Hmatrix_1  , &ione);         // this update is already done right after scf loop 

            memcpy(Kptr[kpt]->Hmatrix_0_cpu, Kptr[kpt]->Hmatrix_1_cpu, matrix_size);

            if(ct.tddft_mode == VECTOR_POT )
            {
                std::complex<double> tem_x = zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->Pxmatrix_cpu, &ione);
                current[0] += std::real(tem_x) * Kptr[kpt]->kp.kweight;
                std::complex<double> tem_y = zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->Pymatrix_cpu, &ione);
                current[1] += std::real(tem_y) * Kptr[kpt]->kp.kweight;
                std::complex<double> tem_z = zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->Pzmatrix_cpu, &ione);
                current[2] += std::real(tem_z) * Kptr[kpt]->kp.kweight;
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, current, 3, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
        Sp->ScalapackBlockAllreduce(current, 3);
        Rmg_Symm->symm_vec(current);

        if(ct.BerryPhase && ct.tddft_mode == VECTOR_POT)
        {
            tot_bp_pol = 0.0;
            for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
                std::complex<double> tem_x = zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->BP_Xml, &ione);
                tot_bp_pol += std::real(tem_x) * Kptr[kpt]->kp.kweight;
            }
            MPI_Allreduce(MPI_IN_PLACE, &tot_bp_pol, 1, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
            Sp->ScalapackBlockAllreduce(&tot_bp_pol, 1);
            //Rmg_BP->CalcBP_tddft(Kptr, tot_bp_pol, matrix_glob, *Sp);
        }

        if(pct.kstart == 0 && pct.gridpe == 0)
        {
            if(ct.tddft_mode == VECTOR_POT )
            {
                fprintf(current_fi, "\n  %f  %18.10e  %18.10e  %18.10e ",
                        tot_steps*time_step, current[0], current[1], current[2]);
                if(ct.BerryPhase) fprintf(dbp_fi, "\n  %f  %18.10e  %18.10e  %18.10e ",
                        tot_steps*time_step, tot_bp_pol, 0.0,0.0);
            }
            else
            {
                fprintf(dfi, "\n  %f  %18.10e  %18.10e  %18.10e ",
                        tot_steps*time_step, dipole_tot[0], dipole_tot[1], dipole_tot[2]);
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
        }

        delete RT2a;

        if((tddft_steps +1) % ct.checkpoint == 0)
        {   
            RT2a = new RmgTimer("2-TDDFT: Write");

            my_sync_device();
            for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
            {
                int kpt_glob = kpt + pct.kstart;

                char newname[MAX_PATH + 20];
                sprintf (newname, "%s_spin%d_kpt%d_gridpe%d", ct.outfile_tddft, pct.spinpe, kpt_glob, pct.gridpe);
                WriteData_rmgtddft(newname, vh.data(), vxc.data(), vh_dipole.data(), (double *)Kptr[kpt]->Pn0_cpu, (double *)Kptr[kpt]->Hmatrix_cpu, 
                        (double *)Kptr[kpt]->Hmatrix_m1_cpu, (double *)Kptr[kpt]->Hmatrix_0_cpu, tot_steps+1, n2, n2_C, Eterms, Hcore_tddft, numst);
            }

            if(pct.kstart == 0 && pct.gridpe == 0)
            {
                if(ct.tddft_mode == VECTOR_POT )
                    fflush(current_fi);
                else
                {
                    fflush(dfi);
                    fflush(efi);
                }
                if(ct.BerryPhase)
                    fflush(dbp_fi);
            }
            delete RT2a;
        }


    }

#if CUDA_ENABLED || HIP_ENABLED
    gpuFree(Kptr[0]->work_dev);
    gpuFree(Kptr[0]->psi_dev);
#endif
    if(pct.kstart == 0 && pct.gridpe == 0)
    {
        if(ct.tddft_mode == VECTOR_POT )
            fclose(current_fi);
        else
        {
            fclose(dfi);
            fclose(efi);
        }
        if(ct.BerryPhase)
            fclose(dbp_fi);
    }

    RT2a = new RmgTimer("2-TDDFT: Write");
    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        int kpt_glob = kpt + pct.kstart;

        char newname[MAX_PATH + 20];
        sprintf (newname, "%s_spin%d_kpt%d_gridpe%d", ct.outfile_tddft, pct.spinpe, kpt_glob, pct.gridpe);
        WriteData_rmgtddft(newname, vh.data(), vxc.data(), vh_dipole.data(), (double *)Kptr[kpt]->Pn0_cpu, (double *)Kptr[kpt]->Hmatrix_cpu, 
                (double *)Kptr[kpt]->Hmatrix_m1_cpu, (double *)Kptr[kpt]->Hmatrix_0_cpu, tot_steps+1, n2, n2_C, Eterms, Hcore_tddft, numst);
    }
    delete RT2a;
    delete RT0;
}

