/*
 *
 * Copyright 2026 The RMG Project Developers. See the COPYRIGHT file 
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

#include "rmg_tddft.h"
#include "../Headers/prototypes_tddft.h"
#include "GatherScatter.h"
#include "rmg_dev_allocate.h"
#include "rmg_hvector.h"
#include "blas_driver.h"

void  init_point_charge_pot(double *vtot_psi, int density);
void eldyn_ort(int *desca, int Mdim, int Ndim, double *F,double *Po0,double *Po1,int *p_Ieldyn,  double *thrs,int*maxiter,  double *errmax,int
        *niter , int *p_iprint, MPI_Comm comm) ;
void eldyn_ort(int *desca, int Mdim, int Ndim, std::complex<double> *F,std::complex<double> *Po0,std::complex<double> *Po1,int *p_Ieldyn,  double *thrs,int*maxiter,  double *errmax,int
        *niter , int *p_iprint, MPI_Comm comm) ;
void eldyn_nonort(int *p_N, double *S, double *F,double *Po0,double *Pn1,int *p_Ieldyn,  double *thrs,int*maxiter,  double *errmax,int *niter , int *p_iprint) ;



template void rmg::tddft<double, double>::tddft_md(void);
template void rmg::tddft<double, std::complex<double>>::tddft_md(void);
template void rmg::tddft<std::complex<double>, std::complex<double>>::tddft_md(void);
template rmg::tddft<double, double>::~tddft(void);
template rmg::tddft<double, std::complex<double>>::~tddft(void);
template rmg::tddft<std::complex<double>, std::complex<double>>::~tddft(void);

template void rmg::tddft<double, double>::gather_rho_matrix(double *, double *);
template void rmg::tddft<double, std::complex<double>>::gather_rho_matrix(double *, std::complex<double> *);
template void rmg::tddft<std::complex<double>, std::complex<double>>::gather_rho_matrix(std::complex<double> *, std::complex<double> *);

template rmg::tddft<double, double>::tddft(spinobj<double> &vxc_in,
             fgobj<double> &vh_in,
             fgobj<double> &vnuc_in,
             spinobj<double> &rho_in,
             fgobj<double> &rhocore_in,
             fgobj<double> &rhoc_in,
             Kpoint<double> **Kptr_in);

template rmg::tddft<double, std::complex<double>>::tddft(spinobj<double> &vxc_in,
             fgobj<double> &vh_in,
             fgobj<double> &vnuc_in,
             spinobj<double> &rho_in,
             fgobj<double> &rhocore_in,
             fgobj<double> &rhoc_in,
             Kpoint<double> **Kptr_in);

template rmg::tddft<std::complex<double>, std::complex<double>>::tddft(spinobj<double> &vxc_in,
             fgobj<double> &vh_in,
             fgobj<double> &vnuc_in,
             spinobj<double> &rho_in,
             fgobj<double> &rhocore_in,
             fgobj<double> &rhoc_in,
             Kpoint<std::complex<double>> **Kptr_in);

template <typename OrbitalType, typename MatrixType>
rmg::tddft<OrbitalType, MatrixType>::tddft(spinobj<double> &vxc_in,
             fgobj<double> &vh_in,
             fgobj<double> &vnuc_in,
             spinobj<double> &rho_in,
             fgobj<double> &rhocore_in,
             fgobj<double> &rhoc_in,
             Kpoint<OrbitalType> **Kptr_in) : vxc(vxc_in), vh(vh_in), vnuc(vnuc_in),
                                            rho(rho_in), rhocore(rhocore_in), rhoc(rhoc_in)
{
    RmgTimer RT0("2-TDDFT: Initialization");
    Kptr = Kptr_in;

#if USE_NCCL
    if(pct.nccl_rank == 0)
    {
        rmg::error(ncclGetUniqueId(&ct.nccl_nd_id));
    }
    MPI_Bcast(&ct.nccl_nd_id, sizeof(ct.nccl_nd_id), MPI_BYTE, 0, pct.nccl_comm);
#if CUDA_ENABLED 
    rmg::error(cuDeviceGet( &ct.cu_dev, 0 ));
    rmg::error(cudaSetDevice(ct.cu_dev));
#endif
    rmg::error(ncclCommInitRank(&ct.nccl_local_comm, pct.nccl_comm_npes, ct.nccl_nd_id, pct.nccl_rank));
#endif  
    iprint = ct.verbose;
    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    pbasis = Kptr[0]->pbasis;
    pbasis_noncoll = pbasis * ct.noncoll_factor;

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
            rmg::error("Invalid subdiag_driver type in TDDFT");

    } // end switch

    if(ct.tddft_gpu)
    {
        scalapack_groups = pct.grid_npes;
    }
    // all tddft propergating will use 1 gpu only, not sure the speed comparison with scalapack for a large system 
    int last = 1;
    numst = ct.num_states - ct.tddft_start_state; 
    MPI_Comm_size(pct.local_comm, &pct.local_comm_npes);
    numst = ( numst/pct.local_comm_npes ) * pct.local_comm_npes; 
    if(ct.tddft_tiledMM == 1)
    {
        if(!ct.tddft_gpu)
        {
            pct.local_comm = pct.grid_comm;
        }
        else
        {
            pct.local_comm = pct.nccl_comm;
        }

        MPI_Comm_size(pct.local_comm, &pct.local_comm_npes);
        MPI_Comm_rank(pct.local_comm, &pct.local_rank);

        eldyn_comm = pct.local_comm;

        // reduce the number of unoccupied states so that numst is divisible by local_comm_npes
        numst = ( numst/pct.local_comm_npes ) * pct.local_comm_npes; 
        n2 = numst * numst/pct.local_comm_npes;
        Mdim = numst;
        Ndim = numst/pct.local_comm_npes;

        Sp = new Scalapack(scalapack_groups, pct.thisimg, ct.images_per_node, numst,
            ct.scalapack_block_factor, last, pct.grid_comm);
    }
    else
    {
        Sp = new Scalapack(scalapack_groups, pct.thisimg, ct.images_per_node, numst,
            ct.scalapack_block_factor, last, pct.grid_comm);
        Mdim = Sp->GetDistMdim();
        Ndim = Sp->GetDistNdim();
        n2 = Sp->GetDistMdim() * Sp->GetDistNdim();

        if(!Sp->Participates()) n2 = 1;
        eldyn_comm =  Sp->GetComm() ;
    }
    int *desca = Sp->GetDistDesca();

    n22 = 2* n2;
    n2_C = n2 * sizeof(MatrixType)/sizeof(double);

    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        Kptr[kpt]->Hmatrix_cpu     = (void *)RmgMallocHost((size_t)n2*sizeof(MatrixType));
        Kptr[kpt]->Pn0_cpu         = (void *)RmgMallocHost((size_t)n2*sizeof(double)*2);
        Kptr[kpt]->Pn1_cpu         = (void *)RmgMallocHost((size_t)n2*sizeof(double)*2);
        Kptr[kpt]->Hmatrix_m1_cpu  = (void *)RmgMallocHost((size_t)n2*sizeof(MatrixType));
        Kptr[kpt]->Hmatrix_1_cpu  = (void *)RmgMallocHost((size_t)n2*sizeof(MatrixType));
        Kptr[kpt]->Hmatrix_0_cpu   = (void *)RmgMallocHost((size_t)n2*sizeof(MatrixType));
        if(ct.tddft_mode == VECTOR_POT)
        {
            Kptr[kpt]->Pxmatrix_cpu   = (std::complex<double> *)RmgMallocHost((size_t)n2*sizeof(std::complex<double>));
            Kptr[kpt]->Pymatrix_cpu   = (std::complex<double> *)RmgMallocHost((size_t)n2*sizeof(std::complex<double>));
            Kptr[kpt]->Pzmatrix_cpu   = (std::complex<double> *)RmgMallocHost((size_t)n2*sizeof(std::complex<double>));

        }
        else
        {
            Kptr[kpt]->Akick_cpu   = (OrbitalType *)RmgMallocHost((size_t)n2*sizeof(OrbitalType));
        }
    }

    matrix_size = n2*sizeof(MatrixType);
#if CUDA_ENABLED || HIP_ENABLED
    rmg_device_pool->malloc(&Hmatrix, n2);
    rmg_device_pool->malloc(&Hmatrix_m1, n2);
    rmg_device_pool->malloc(&Hmatrix_0, n2);
    rmg_device_pool->malloc(&Hmatrix_1, n2);
    if(typeid(MatrixType) == typeid(double) )
    {
        rmg_device_pool->malloc(&Pn0, 2*n2);
        rmg_device_pool->malloc(&Pn1, 2*n2);
    }
    else
    {
        rmg_device_pool->malloc(&Pn0, n2);
        rmg_device_pool->malloc(&Pn1, n2);
    }


    size_t num_psik = ct.num_kpts_pe * ct.num_states * (size_t)pbasis_noncoll;
    rmg_device_pool->malloc(&psi_dev_pool, num_psik * 2);
    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        Kptr[kpt]->psi_dev  = psi_dev_pool + 2 * kpt * ct.num_states * pbasis_noncoll;
        Kptr[kpt]->work_dev = Kptr[kpt]->psi_dev + ct.num_states * pbasis_noncoll;
        // even for tddft_floatprecision, psi_dev is still used in VecPHmatrix and CurrentNlpp
        RmgMemcpy(Kptr[kpt]->psi_dev, Kptr[kpt]->orbital_storage, ct.num_states * pbasis_noncoll * sizeof(OrbitalType));

        if(ct.tddft_floatprecision)
        {

            size_t psi_alloc = (size_t)ct.num_states * (size_t)pbasis_noncoll * sizeof(OrbitalType);
            gpuMalloc((void **)&Kptr[kpt]->psi_dev_float, psi_alloc/2);
            gpuMalloc((void **)&Kptr[kpt]->work_dev_float, psi_alloc/2);

            size_t count = (size_t)ct.num_states * (size_t)pbasis_noncoll;
            if(typeid(OrbitalType) == typeid(double))
            {
                float *work_conv = new float[count];
                CopyAndConvert(count, (double *)Kptr[kpt]->orbital_storage, work_conv);
                RmgMemcpy(Kptr[kpt]->psi_dev_float, work_conv, count * sizeof(float));
                delete [] work_conv;
            }
            else if(typeid(OrbitalType) == typeid(std::complex<double>))
            {
                std::complex<float> *work_conv = new std::complex<float>[count];
                CopyAndConvert(count, (std::complex<double> *)Kptr[kpt]->orbital_storage, work_conv);
                RmgMemcpy(Kptr[kpt]->psi_dev_float, work_conv, count * sizeof(std::complex<float>));
                delete [] work_conv;
            }

        }
    }

#else
    Pn1  = (MatrixType *)RmgMallocHost((size_t)n2*sizeof(double)*2);
    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        Kptr[kpt]->work_cpu = new OrbitalType[(size_t)ct.num_states * (size_t)pbasis_noncoll];
    }
#endif

    std::vector<double> diag_elem(numst);

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

    efactor = ct.energy_output_conversion[ct.energy_output_units];
    eunits = ct.energy_output_string[ct.energy_output_units].c_str();

    if(!ct.norm_conserving_pp)
    {
        rmg::error(" \n  TDDFT support NCPP only \n");
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
        }
    }


    double time_step = ct.tddft_time_step;

    ReadData (ct.infile, vh.data(), rho_ground.data(), vxc.data(), Kptr);
    rho_ground.get_oppo();

    if(ct.tddft_energy)
    {

        tddft_energy_init(vxc, vh, vnuc, rho_ground, rhocore, rhoc, Kptr, *Sp, Mdim, Ndim, Eterms_ground);
        for(int i = 0; i < 6; i++) Eterms[i] = Eterms_ground[i];
        if(pct.kstart == 0 && pct.gridpe == 0 && pct.spinpe == 0)
        {
            filename = std::string(ct.basename) + "_totalE";
            efi = fopen(filename.c_str(), "w");

            fprintf(efi, " && totalE_0, EkinPseudo_0, Vh_0, Exc_0  II  Edownfold%s", eunits);
            fprintf(efi, "\n&& %16.8e  %16.8e  %16.8e  %16.8e %16.8e %16.8e at Ground state", Eterms_ground[0], 
                    Eterms_ground[1],Eterms_ground[2],Eterms_ground[3],Eterms_ground[4],Eterms_ground[5]);
        }
    }

    if(ct.restart_tddft)
    {

        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {
            int kpt_glob = kpt + pct.kstart;

            std::string ofile = std::format("{}_spin{}_kpt{}_gridpe{}", ct.infile_tddft, pct.spinpe, kpt_glob, pct.gridpe);
            ReadData_rmgtddft(ofile.c_str(), vh.data(), vxc.data(), vh_dipole.data(), (double *)Kptr[kpt]->Pn0_cpu, (double *)Kptr[kpt]->Hmatrix_cpu, 
                    (double *)Kptr[kpt]->Hmatrix_m1_cpu, (double *)Kptr[kpt]->Hmatrix_0_cpu, 
                    &pre_steps, n2, n2_C, numst);
        }
    }
    else
    {

        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {

            for(int i = 0; i < numst; i++) diag_elem[i] = Kptr[kpt]->Kstates[i + ct.tddft_start_state].eig[0];

            double one = 1.0;
            memset(Kptr[kpt]->Hmatrix_cpu, 0, n2*sizeof(MatrixType));
            MatDiagSet((MatrixType *)Kptr[kpt]->Hmatrix_cpu, diag_elem, one, numst, *Sp);

            memcpy(Kptr[kpt]->Hmatrix_1_cpu, Kptr[kpt]->Hmatrix_cpu, matrix_size);
            memcpy(Kptr[kpt]->Hmatrix_0_cpu, Kptr[kpt]->Hmatrix_cpu, matrix_size);
            memcpy(Kptr[kpt]->Hmatrix_m1_cpu, Kptr[kpt]->Hmatrix_0_cpu, matrix_size);


            for(int i = 0; i < numst; i++) diag_elem[i] =  Kptr[kpt]->Kstates[i + ct.tddft_start_state].occupation[0];
            memset(Kptr[kpt]->Pn0_cpu, 0, 2*n2*sizeof(double));
            MatDiagSet((MatrixType *)Kptr[kpt]->Pn0_cpu, diag_elem, one, numst, *Sp);
        }

    }

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
            HmatrixUpdate<OrbitalType, OrbitalType, OrbitalType>(Kptr[kpt], vtot_psi, vxc_psi, (OrbitalType *)Kptr[kpt]->Akick_cpu, ct.tddft_start_state, numst, desca);
            daxpy ( &n2_C,  &alpha, (double *)Kptr[kpt]->Akick_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_0_cpu , &ione) ;
            memcpy(Kptr[kpt]->Hmatrix_m1_cpu, Kptr[kpt]->Hmatrix_0_cpu, matrix_size);
        }
    }

    //  initialize   data for rt-td-dft
    //int nblock = 10 ;   //  size of tthe block for printing (debug!)

    /*
       if(pct.gridpe == 0) { printf("**** Hmat  : \n");  print_matrix_d(Hmatrix,   &nblock, &numst)   ; }
       if(pct.gridpe == 0) { printf("**** Hmat0 : \n");  print_matrix_d(Hmatrix_0, &nblock, &numst)   ; }
       if(pct.gridpe == 0) { printf("**** Hmat1 : \n");  print_matrix_d(Hmatrix_1, &nblock, &numst)   ; }
     */

    current[0] = 0.0;
    current[1] = 0.0;
    current[2] = 0.0;
    current0[0] = 0.0;
    current0[1] = 0.0;
    current0[2] = 0.0;
    tot_bp_pol = 0.0;
    if(ct.tddft_mode == VECTOR_POT)
    {
        // vector potential will be A(t) =  ct.efield_tddft * cos(tddft_frequency * t)
        // VecP matrix is <psi| ct.efied_tddft dot gradient | psi> 
        //
        if(ct.verbose) {
            rmg::printlog("\n starting VecP matrix ");
            fflush(NULL);
        }
        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
            VecPHmatrix(Kptr[kpt], ct.efield_tddft_crds, desca, ct.tddft_start_state, numst);

            if(pre_steps == 0)
            {
                // at t= 0, cos(omega t) = 1.0
                daxpy ( &n2_C ,  &ct.efield_tddft_crds[0], (double *)Kptr[kpt]->Pxmatrix_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_m1_cpu,  &ione) ;
                daxpy ( &n2_C ,  &ct.efield_tddft_crds[1], (double *)Kptr[kpt]->Pymatrix_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_m1_cpu,  &ione) ;
                daxpy ( &n2_C ,  &ct.efield_tddft_crds[2], (double *)Kptr[kpt]->Pzmatrix_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_m1_cpu,  &ione) ;
                memcpy(Kptr[kpt]->Hmatrix_0_cpu, Kptr[kpt]->Hmatrix_m1_cpu, matrix_size);
            }

            CurrentNlpp(Kptr[kpt], desca, ct.tddft_start_state, numst);

            if(0)
            {
                for(int i = 0; i < n2; i++) 
                {
                    Kptr[kpt]->Pxmatrix_cpu[i] = 0.0;
                    Kptr[kpt]->Pymatrix_cpu[i] = 0.0;
                    Kptr[kpt]->Pzmatrix_cpu[i] = 0.0;
                }
                CurrentOperator(Kptr[kpt], desca, ct.tddft_start_state);
            }
        }
        if(ct.verbose) {
            rmg::printlog("\n done VecP matrix ");
            fflush(NULL);
        }

        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
            std::complex<double> tem_x = rmg_zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->Pxmatrix_cpu, &ione);
            current0[0] += std::real(tem_x) * Kptr[kpt]->kp.kweight;
            std::complex<double> tem_y = rmg_zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->Pymatrix_cpu, &ione);
            current0[1] += std::real(tem_y) * Kptr[kpt]->kp.kweight;
            std::complex<double> tem_z = rmg_zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->Pzmatrix_cpu, &ione);
            current0[2] += std::real(tem_z) * Kptr[kpt]->kp.kweight;
        }
        if(ct.BerryPhase)
        {
            // Rmg_BP->CalcBP_Skk1(Kptr, ct.tddft_start_state, matrix_glob, *Sp);
            // Rmg_BP->CalcBP_tddft(Kptr, tot_bp_pol, matrix_glob, *Sp);
            Rmg_BP->tddft_Xml(Kptr, ct.tddft_start_state, *Sp);
            tot_bp_pol = 0.0;
            for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
                std::complex<double> tem_x = rmg_zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->BP_Xml, &ione);
                tot_bp_pol += std::real(tem_x) * Kptr[kpt]->kp.kweight;
            }
            MPI_Allreduce(MPI_IN_PLACE, &tot_bp_pol, 1, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
            MPI_Allreduce(MPI_IN_PLACE, &tot_bp_pol, 1, MPI_DOUBLE, MPI_SUM, eldyn_comm);
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, current0, 3, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
    MPI_Allreduce(MPI_IN_PLACE, current0, 3, MPI_DOUBLE, MPI_SUM, eldyn_comm);
    Rmg_Symm->symm_vec(current0);

#if CUDA_ENABLED || HIP_ENABLED
    if(ct.tddft_floatprecision)
    {
        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {

            // psi_dev is no longer needed, so psi_dev_float and work_dev_float use these memory
            Kptr[kpt]->psi_dev_float = Kptr[kpt]->psi_dev;
            Kptr[kpt]->work_dev_float = Kptr[kpt]->work_dev;

            size_t count = (size_t)ct.num_states * (size_t)pbasis_noncoll;
            if(typeid(OrbitalType) == typeid(double))
            {
                float *work_conv = new float[count];
                CopyAndConvert(count, (double *)Kptr[kpt]->orbital_storage, work_conv);
                RmgMemcpy(Kptr[kpt]->psi_dev_float, work_conv, count * sizeof(float));
                delete [] work_conv;
            }
            else if(typeid(OrbitalType) == typeid(std::complex<double>))
            {
                std::complex<float> *work_conv = new std::complex<float>[count];
                CopyAndConvert(count, (std::complex<double> *)Kptr[kpt]->orbital_storage, work_conv);
                RmgMemcpy(Kptr[kpt]->psi_dev_float, work_conv, count * sizeof(std::complex<float>));
                delete [] work_conv;
            }

        }
    }
#endif


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
        fflush(NULL);
    }

}

    // TDDFT MD loop
template <typename OrbitalType, typename MatrixType>
void rmg::tddft<OrbitalType, MatrixType>::tddft_md(void)
{

    RmgTimer *RT2a;
    Kpoint<double> *kptr_d;
    Kpoint<std::complex<double>> *kptr_c;
    int ij_err=0;
    double vtxc, etxc;
    int *desca = Sp->GetDistDesca();
    static double total_time = 0.0;

    // Recompute sint arrays which is necessary for dynamics.
    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
#if HIP_ENABLED || CUDA_ENABLED
        this->Kptr[kpt]->BetaProjector->project(Kptr[kpt], Kptr[kpt]->newsint_local, 0,
                Kptr[kpt]->nstates * ct.noncoll_factor, Kptr[kpt]->nl_weight_gpu);
#else
        this->Kptr[kpt]->BetaProjector->project(Kptr[kpt], Kptr[kpt]->newsint_local, 0,
                Kptr[kpt]->nstates * ct.noncoll_factor, Kptr[kpt]->nl_weight);
#endif
    }

    //  run rt-td-dft
    for(int tddft_steps = 0; tddft_steps < ct.tddft_steps; tddft_steps++)
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

        rmg::sync_device();
        delete RT2a;


        int  Max_iter_scf = 10 ; int  iter_scf =0 ;
        double err =1.0e0   ;
        thrs_dHmat  = 1e-7  ;

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
                    Hmatrix = (MatrixType *)Kptr[kpt]->Hmatrix_cpu;
                    Hmatrix_m1 = (MatrixType *)Kptr[kpt]->Hmatrix_m1_cpu;
                    Hmatrix_1 = (MatrixType *)Kptr[kpt]->Hmatrix_1_cpu;
                    Hmatrix_0 = (MatrixType *)Kptr[kpt]->Hmatrix_0_cpu;
                    Pn0 = (MatrixType *)Kptr[kpt]->Pn0_cpu;
                    Pn1 = (MatrixType *)Kptr[kpt]->Pn1_cpu;
                }
                rmg::sync_device();
                delete RT2a;
                RT2a = new RmgTimer("2-TDDFT: ELDYN");
                if(ct.verbose) {
                    rmg::printlog("\n start magnus and eldyn ");
                    fflush(NULL);
                }
                magnus ((double *)Hmatrix_0,    (double *)Hmatrix_1 , time_step, (double *)Hmatrix_m1 , n2_C) ; 
                /* --- C++  version:  --*/

                eldyn_ort(desca, Mdim, Ndim,  Hmatrix_m1,Pn0,Pn1,&Ieldyn, &thrs_bch,&maxiter_bch,  &errmax_bch,&niter_bch ,  &iprint, eldyn_comm) ;
                RmgMemcpy(Kptr[kpt]->Pn1_cpu, Pn1, 2*n2*sizeof(double));

                if(ct.verbose) {
                    rmg::printlog("\n done magnus and eldyn ");
                    fflush(NULL);
                }
                delete(RT2a);

                // if(pct.gridpe == 0) { printf("**** Pn1 : \n");   print_matrix_z(Pn1,  &nblock, &numst)  ; }

                //            for(i = 0; i < 10; i++) 
                //            { printf("Pn\n");
                //           for(int j = 0; j < 10; j++) printf(" %8.1e", i, Pn1[i*numst + j]);
                //          }


                /////// <----- update Hamiltonian from  Pn1
                RT2a = new RmgTimer("2-TDDFT: Rho");
                rmg::sync_device();
                if(ct.verbose) {
                    rmg::printlog("\n start rho calc ");
                    fflush(NULL);
                }
                if(ct.tddft_floatprecision)
                {
                    if(ct.is_gamma)
                    {
                        kptr_d = (Kpoint<double> *)Kptr[kpt];
                        GetNewRho_rmgtddft<double, float, MatrixType>(kptr_d, rho_k, Pn1, numst, ct.tddft_start_state, *Sp);
                    }
                    else
                    {
                        kptr_c = (Kpoint<std::complex<double>> *)Kptr[kpt];
                        GetNewRho_rmgtddft<std::complex<double>, std::complex<float>, std::complex<double> >(kptr_c, rho_k, (std::complex<double> *)Pn1, numst, ct.tddft_start_state, *Sp);
                    }
                }
                else
                {
                    GetNewRho_rmgtddft<OrbitalType, OrbitalType, MatrixType>(Kptr[kpt], rho_k, Pn1, numst, ct.tddft_start_state, *Sp);
                }

                if(ct.verbose) {
                    rmg::printlog("\n done rho calc ");
                    fflush(NULL);
                }
                int kpt_glob = kpt + pct.kstart;
                for(int idx = 0; idx < FP0_BASIS; idx++) rho_ksum[idx] += rho_k[idx] * ct.kp[kpt_glob].kweight;

                delete(RT2a);

            }

            MPI_Allreduce(MPI_IN_PLACE, rho_ksum.data(), FP0_BASIS, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
            for(int idx = 0; idx < FP0_BASIS; idx++) rho[idx] = rho_ksum[idx] + rho_ground[idx];
            rho.get_oppo();


            //write_rho_x(rho, "update rho");
            //write_rho_x(rho_ground.data(), "groumd rho");

            //dcopy(&FP0_BASIS, vh_dipole.data(), &ione, vh_dipole_old.data(), &ione);
            //dcopy(&FP0_BASIS, vh.data(), &ione, vh_old.data(), &ione);
            //dcopy(&FP0_BASIS, vxc.data(), &ione, vxc_old.data(), &ione);
            vh_old = vh;
            vh_dipole_old = vh_dipole;
            vxc_old = vxc;

            //get_vxc(rho, rho_oppo, rhocore, vxc);
            RmgTimer *RT1 = new RmgTimer("2-TDDFT: exchange/correlation");
            //vxc_in = vxc;
            compute_vxc(rho.data(), rhocore.data(), etxc, vtxc, vxc.data(), ct.nspin);
            delete RT1;

            RT1 = new RmgTimer("2-TDDFT: Vh");
            //vh_in = vh;
            VhDriver(rho.data(), rhoc.data(), vh.data(), ct.vh_ext, 1.0-12);
            delete RT1;

            get_dipole(rho.data(), rhoc.data(), dipole_tot);
            if(ct.dipole_corr[0]+ct.dipole_corr[1]+ct.dipole_corr[2] >0)
            {
                DipoleCorrection(dipole_tot,  vh_dipole.data());
            }

            // noncoll need change
            for (int idx = 0; idx < FP0_BASIS; idx++) {
                vtot[idx] = vxc[idx] + vh[idx] + vh_dipole[idx]
                    -vxc_old[idx] -vh_old[idx] - vh_dipole_old[idx];
            }

            vxc_diff = vxc - vxc_old;


            //noncoll need change
            GetVtotPsi (vtot_psi.data(), vtot.data(), Rmg_G->default_FG_RATIO);

            for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
                RT2a = new RmgTimer("2-TDDFT: Hupdate");
                if(ct.tddft_floatprecision)
                {
                    if(ct.is_gamma)
                    {
                        kptr_d = (Kpoint<double> *)Kptr[kpt];
                        HmatrixUpdate<double,float, MatrixType>  (kptr_d, vtot_psi, vxc_psi, (MatrixType *)Kptr[kpt]->Hmatrix_m1_cpu, ct.tddft_start_state, numst, desca);                                     
                    }
                    else
                    {
                        kptr_c = (Kpoint<std::complex<double>> *)Kptr[kpt];
                        HmatrixUpdate<std::complex<double>,std::complex<float>, std::complex<double>>  (kptr_c, vtot_psi, vxc_psi, (std::complex<double> *)Kptr[kpt]->Hmatrix_m1_cpu, ct.tddft_start_state, numst, desca);                                     
                    }
                }
                else
                {
                    HmatrixUpdate<OrbitalType,OrbitalType, MatrixType>  (Kptr[kpt], vtot_psi, vxc_psi, (MatrixType *)Kptr[kpt]->Hmatrix_m1_cpu, ct.tddft_start_state, numst, desca);                                     
                }
                delete(RT2a);

                RT2a = new RmgTimer("2-TDDFT: conv check");
                rmg::sync_device();
                double one = 1.0, mone = -1.0;
                daxpy( &n2_C ,  &one, (double *)Kptr[kpt]->Hmatrix_m1_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_cpu,  &ione) ;

                //////////  < ---  end of Hamiltonian update

                // check error and update Hmatrix_1:
                rmg::sync_device();
                daxpy ( &n2_C ,  &mone, (double *)Kptr[kpt]->Hmatrix_cpu, &ione , (double *)Kptr[kpt]->Hmatrix_1_cpu ,  &ione) ;

                //tst_conv_matrix (&err, &ij_err ,  Hmatrix_1,  n2, Sp->GetComm()) ;  //  check error  how close  H and H_old are

                bool tConv;
                tstconv((double *)Kptr[kpt]->Hmatrix_1_cpu, &n2_C, &thrs_dHmat,&ij_err,&err,&tConv, eldyn_comm);
                memcpy(Kptr[kpt]->Hmatrix_1_cpu, Kptr[kpt]->Hmatrix_cpu, matrix_size);
                delete(RT2a);
            }

            MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, pct.kpsub_comm);
            MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, pct.spin_comm);


            if(pct.imgpe == 0) { printf("step: %5d  iteration: %d  thrs= %12.5e err=  %12.5e at element: %5d \n", 
                    tddft_steps, iter_scf,    thrs_dHmat,  err,         ij_err); } 
            rmg::printlog("step: %5d  iteration: %d  thrs= %12.5e err=  %12.5e at element: %5d \n", 
                    tddft_steps, iter_scf,    thrs_dHmat,  err,         ij_err);  
            //err= -1.0e0 ;  
            iter_scf ++ ;
        } //---- end of  SCF/while loop 


        RT2a = new RmgTimer("2-TDDFT: current and dipole");
        //  extract dipole from rho(Pn1)
        get_dipole(rho.data(), rhoc.data(), dipole_tot);
        /*  done with propagation,  save Pn1 ->  Pn0 */
        if(ct.tddft_energy)
        {
            static int header_once;
            Eterms[3] = etxc;
            tddft_energy(vh, rho, rhoc, Kptr, Mdim, Ndim, Eterms, eldyn_comm);
            if(tot_steps == 0)
            {
                if(pct.kstart == 0 && pct.gridpe == 0 && pct.spinpe == 0 && header_once == 0)
                {
                    Eterms_1step = Eterms;
                    fprintf(efi, "\n&& %16.8e  %16.8e  %16.8e  %16.8e %16.8e %16.8e at 1st TDDFT step", Eterms_1step[0],
                            Eterms_1step[1],Eterms_1step[2],Eterms_1step[3],Eterms_1step[4],Eterms_1step[5]);
                    fprintf(efi, "\n&&time  totalE_diff, EkinPseudo_diff, Vh_diff, Exc_diff  %s", eunits);
                    header_once++;
                }
                //for(int i = 0; i < 6; i++) Eterms_1step[i] = Eterms[i];
            }
            if(pct.kstart == 0 && pct.gridpe == 0 && pct.spinpe == 0)
            {
                fprintf(efi, "\n  %f  %16.8e %16.8e,%16.8e,%16.8e   ",
                        total_time, (Eterms[0] - Eterms_1step[0]) * efactor, (Eterms[1] - Eterms_1step[1]) * efactor, 
                        (Eterms[2] - Eterms_1step[2]) * efactor, (Eterms[3] - Eterms_1step[3]) * efactor);
            }
        }


        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
            memcpy(Kptr[kpt]->Pn0_cpu, Kptr[kpt]->Pn1_cpu, n22 * sizeof(double));

            // save current  H0, H1 for the  next step extrapolatiion
            memcpy(Kptr[kpt]->Hmatrix_m1_cpu, Kptr[kpt]->Hmatrix_0_cpu, matrix_size);
            //dcopy(&n2, Hmatrix  , &ione, Hmatrix_1  , &ione);         // this update is already done right after scf loop 

            memcpy(Kptr[kpt]->Hmatrix_0_cpu, Kptr[kpt]->Hmatrix_1_cpu, matrix_size);

            if(ct.tddft_mode == VECTOR_POT )
            {
                std::complex<double> tem_x = rmg_zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->Pxmatrix_cpu, &ione);
                current[0] += std::real(tem_x) * Kptr[kpt]->kp.kweight;
                std::complex<double> tem_y = rmg_zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->Pymatrix_cpu, &ione);
                current[1] += std::real(tem_y) * Kptr[kpt]->kp.kweight;
                std::complex<double> tem_z = rmg_zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->Pzmatrix_cpu, &ione);
                current[2] += std::real(tem_z) * Kptr[kpt]->kp.kweight;
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, current, 3, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
        MPI_Allreduce(MPI_IN_PLACE, current, 3, MPI_DOUBLE, MPI_SUM, eldyn_comm);
        Rmg_Symm->symm_vec(current);

        if(ct.BerryPhase && ct.tddft_mode == VECTOR_POT)
        {
            tot_bp_pol = 0.0;
            for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
                std::complex<double> tem_x = rmg_zdotc(&n2, (std::complex<double> *)Kptr[kpt]->Pn0_cpu, &ione, (std::complex<double> *)Kptr[kpt]->BP_Xml, &ione);
                tot_bp_pol += std::real(tem_x) * Kptr[kpt]->kp.kweight;
            }
            MPI_Allreduce(MPI_IN_PLACE, &tot_bp_pol, 1, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
            MPI_Allreduce(MPI_IN_PLACE, &tot_bp_pol, 1, MPI_DOUBLE, MPI_SUM, eldyn_comm);
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
            }
        }

        delete RT2a;

        if((tddft_steps +1) % ct.checkpoint == 0)
        {   
            RT2a = new RmgTimer("2-TDDFT: Write");

            rmg::sync_device();
            for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
            {
                int kpt_glob = kpt + pct.kstart;

                std::string ofile = std::format("{}_spin{}_kpt{}_gridpe{}", 
                        ct.outfile_tddft, pct.spinpe, kpt_glob, pct.gridpe);
                WriteData_rmgtddft(ofile.c_str(), vh.data(), vxc.data(), vh_dipole.data(), (double *)Kptr[kpt]->Pn0_cpu, (double *)Kptr[kpt]->Hmatrix_cpu, 
                        (double *)Kptr[kpt]->Hmatrix_m1_cpu, (double *)Kptr[kpt]->Hmatrix_0_cpu, tot_steps+1, n2, n2_C, numst);
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

        total_time += time_step;
    } // end tddft md loop

    /*When running MD, force pointers need to be rotated before calculating new forces */
    if(ct.forceflag == TDDFT_CVE)
    {
        ct.fpt[0] = 0;
        ct.fpt[1] = 1;
        ct.fpt[2] = 2;
        ct.fpt[3] = 3;
        ct.sqrt_interpolation = false;

        spinobj<double> trho;
        trho = rho;

        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            Atoms[ion].RotateForces();
        }
        rmg::hvector<OrbitalType> rho_matrix_global(ct.num_states*ct.num_states); 
        gather_rho_matrix(rho_matrix_global.data(), Pn1);
        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {
            Kptr[kpt]->save_sint();
            rmg::rotate_sint(Kptr[kpt], Kptr[kpt]->newsint_local, rho_matrix_global.data());
        }
        Force (trho.up.data(), trho.dw.data(), rhoc.data(), vh.data(), vh.data(), vxc.data(), vxc.data(), vnuc.data(), Kptr);
        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) Kptr[kpt]->restore_sint();
        get_ddd (vtot.data(), vxc.data(), true);

    }

}

template <typename OrbitalType, typename MatrixType>
rmg::tddft<OrbitalType, MatrixType>::~tddft(void)
{
    RmgTimer *RT2a;
#if CUDA_ENABLED || HIP_ENABLED
    rmg_device_pool->free(psi_dev_pool);
    rmg_device_pool->free(Pn1);
    rmg_device_pool->free(Pn0);
    rmg_device_pool->free(Hmatrix_1);
    rmg_device_pool->free(Hmatrix_0);
    rmg_device_pool->free(Hmatrix_m1);
    rmg_device_pool->free(Hmatrix);
#endif  
    if(pct.kstart == 0 && pct.gridpe == 0)
    {
        if(ct.tddft_mode == VECTOR_POT )
            fclose(current_fi);
        else
        {
            fclose(dfi);
        }
        if(ct.BerryPhase)
            fclose(dbp_fi);
    }
    if(ct.tddft_energy && pct.kstart == 0 && pct.gridpe == 0 && pct.spinpe == 0)
    {
        fclose(efi);
    }

    RT2a = new RmgTimer("2-TDDFT: Write");
    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        int kpt_glob = kpt + pct.kstart;

        std::string ofile = std::format("{}_spin{}_kpt{}_gridpe{}",
                ct.outfile_tddft, pct.spinpe, kpt_glob, pct.gridpe);
        WriteData_rmgtddft(ofile.c_str(), vh.data(), vxc.data(), vh_dipole.data(), (double *)Kptr[kpt]->Pn0_cpu, (double *)Kptr[kpt]->Hmatrix_cpu,
                (double *)Kptr[kpt]->Hmatrix_m1_cpu, (double *)Kptr[kpt]->Hmatrix_0_cpu, tot_steps+1, n2, n2_C, numst);
    }
    delete RT2a;

    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        RmgFreeHost(Kptr[kpt]->Hmatrix_cpu);
        RmgFreeHost(Kptr[kpt]->Pn0_cpu);
        RmgFreeHost(Kptr[kpt]->Pn1_cpu);
        RmgFreeHost(Kptr[kpt]->Hmatrix_m1_cpu);
        RmgFreeHost(Kptr[kpt]->Hmatrix_1_cpu);
        RmgFreeHost(Kptr[kpt]->Hmatrix_0_cpu);
        if(ct.tddft_mode == VECTOR_POT)
        {
            RmgFreeHost(Kptr[kpt]->Pxmatrix_cpu);
            RmgFreeHost(Kptr[kpt]->Pymatrix_cpu);
            RmgFreeHost(Kptr[kpt]->Pzmatrix_cpu);

        }
        else
        {
            RmgFreeHost(Kptr[kpt]->Akick_cpu);
        }
    }

    if(Sp) delete Sp;
}


template <typename OrbitalType, typename MatrixType>
void rmg::tddft<OrbitalType, MatrixType>::tstconv(double *C,int *p_M, double *p_thrs,int *p_ierr, double *p_err, bool *p_tconv, MPI_Comm comm) 
{
    int     M     = *(p_M)    ;  //  [in]  :  total  size of matrix (2*Nbasis*Nbasis)
    double  thrs  = *(p_thrs) ;  //  [in]  :  convergence threshold
    double  err               ;  //  [out] :   error= abs of max element in the matrix
    int    ierr   =   0       ;  //  [out] :   location of err in matrix/vector
    bool   tconv  =  false    ;  //  [out] :   if converged ?  true or false?

    rmg::sync_device();
    if(!ct.tddft_gpu)
    { 
        err = abs(C[0]); 
        for (int i=0; i <M ;i++) {
            double err_tmp = abs(C[i]) ; 
            if (err_tmp > err) {
                err  = err_tmp ;
                ierr = i       ;
            }
        }
    }
    else
    {
#if CUDA_ENABLED || HIP_ENABLED
        int idx;
#if HIP_ENABLED
        hipblasIdamax(ct.gpublas_handle, M, C, 1, &idx);
#endif
#if CUDA_ENABLED
        cublasIdamax(ct.gpublas_handle, M, C, 1, &idx);
#endif
        idx -=1;    
        // hipblasIdamax return the index in fortran way, starting from 1
        gpuMemcpy(&err, &C[idx], sizeof(double), gpuMemcpyDeviceToHost);
        err = abs(err);
        ierr = idx;
#endif
    }

    rmg::sync_device();
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, comm);

    if (err < thrs)  tconv = true ;
    /*-- return values **/
    *(p_err )  =  err  ;
    *(p_ierr)  =  ierr ;
    *(p_tconv) =  tconv ;
}


template <typename OrbitalType, typename MatrixType>
void rmg::tddft<OrbitalType, MatrixType>::tstconv(float *C,int *p_M, double *p_thrs,int *p_ierr, double *p_err, bool *p_tconv, MPI_Comm comm) 
{
    int     M     = *(p_M)    ;  //  [in]  :  total  size of matrix (2*Nbasis*Nbasis)
    double  thrs  = *(p_thrs) ;  //  [in]  :  convergence threshold
    float  err               ;  //  [out] :   error= abs of max element in the matrix
    int    ierr   =   0       ;  //  [out] :   location of err in matrix/vector
    bool   tconv  =  false    ;  //  [out] :   if converged ?  true or false?


    rmg::sync_device();
#if CUDA_ENABLED || HIP_ENABLED
    int idx;
#if HIP_ENABLED
    hipblasIsamax(ct.gpublas_handle, M, C, 1, &idx);
#endif
#if CUDA_ENABLED
    cublasIsamax(ct.gpublas_handle, M, C, 1, &idx);
#endif
    idx -=1;    
    gpuMemcpy(&err, &C[idx], sizeof(float), gpuMemcpyDeviceToHost);
    err = abs(err);
    ierr = idx;
#else
    err = abs(C[0]); 
    for (int i=0; i <M ;i++) {
        double err_tmp = abs(C[i]) ; 
        if (err_tmp > err) {
            err  = err_tmp ;
            ierr = i       ;
        }
    }
#endif

    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_FLOAT, MPI_MAX, comm);

    if (err < thrs)  tconv = true ;
    /*-- return values **/
    *(p_err )  =  (double)err  ;
    *(p_ierr)  =  ierr ;
    *(p_tconv) =  tconv ;
}


template <typename OrbitalType, typename MatrixType>
void rmg::tddft<OrbitalType, MatrixType>::gather_rho_matrix(OrbitalType *rho_matrix_global, MatrixType *rho_matrix)
{   

    if(ct.tddft_tiledMM == 1)
    {
        double *rho_R = (double *)rho_matrix_global;
        std::complex<double> *rho_C = (std::complex<double> *)rho_matrix_global;
        size_t recvcount = numst * numst/pct.local_comm_npes * sizeof(OrbitalType)/sizeof(double);
        if(typeid(OrbitalType) == typeid(double))
        {
            for(int i = 0; i < numst * numst/pct.local_comm_npes; i++)
            {
                rho_R[pct.local_rank * numst * numst/pct.local_comm_npes + i] = std::real(rho_matrix[i]);
            }
        }
        else
        {
            for(int i = 0; i < numst * numst/pct.local_comm_npes; i++)
            {
                rho_C[pct.local_rank * numst * numst/pct.local_comm_npes + i] = rho_matrix[i];
            }
        }
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 
                rho_matrix_global, recvcount, MPI_DOUBLE, pct.local_comm);       
    }
    else
    {
        int Mdim = Sp->GetDistMdim();
        int Ndim = Sp->GetDistNdim();
        std::vector<OrbitalType> rho_matrix_dist(Mdim * Ndim);
        double *rho_R = (double *)rho_matrix_dist.data();
        std::complex<double> *rho_C = (std::complex<double> *)rho_matrix_dist.data();
        if(typeid(OrbitalType) == typeid(double))
        {
            for(int i = 0; i < Mdim * Ndim; i++)
            {
                rho_R[i] = std::real(rho_matrix[i]);
            }
        }
        else
        {
            for(int i = 0; i < Mdim * Ndim; i++)
            {
                rho_C[i] = rho_matrix[i];
            }
        }

        for(int i = 0; i < numst * numst; i++) rho_matrix_global[i] = 0.0;
        Sp->GatherEigvectors(rho_matrix_global, rho_matrix_dist.data());
    }

}
