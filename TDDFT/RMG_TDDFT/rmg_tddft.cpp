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

template rmg::tddft<double, double>::tddft(spinobj<double> &vxc_in,
             fgobj<double> &vh_in,
             fgobj<double> &vnuc_in,
             spinobj<double> &rho_in,
             fgobj<double> &rhocore_in,
             fgobj<double> &rhoc_in,
             Kpoint<double> **Kptr_in);


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

    Kptr = Kptr_in;
#if USE_NCCL
    int nlocal_ranks;  
    MPI_Comm_size(pct.local_comm, &nlocal_ranks);
    if(pct.local_rank == 0)
    {
        rmg::error(ncclGetUniqueId(&ct.nccl_nd_id));
    }
    MPI_Bcast(&ct.nccl_nd_id, sizeof(ct.nccl_nd_id), MPI_BYTE, 0, pct.local_comm);
#if CUDA_ENABLED  
    rmg::error(cuDeviceGet( &ct.cu_dev, 0 ));
    rmg::error(cudaSetDevice(ct.cu_dev));
#endif
    rmg::error(ncclCommInitRank(&ct.nccl_local_comm, nlocal_ranks, ct.nccl_nd_id, pct.local_rank));
#endif

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
    int Mdim, Ndim;
    MPI_Comm_size(pct.local_comm, &pct.local_comm_npes);
    numst = ( numst/pct.local_comm_npes ) * pct.local_comm_npes;
    if(ct.tddft_tiledMM == 1)
    {
        if(!ct.tddft_gpu)
        {
            pct.local_rank = pct.gridpe;
            pct.local_comm = pct.grid_comm;
        }
        MPI_Comm_size(pct.local_comm, &pct.local_comm_npes);

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
    desca = Sp->GetDistDesca();
    n22 = 2* n2;
    int n2_C = n2 * sizeof(MatrixType)/sizeof(double);
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

    std::vector<double> diag_elem(numst);
    // Jacek: 
    //double *dHmatrix    = new double[n2];   // storage for  H1 -H1_old 
    MatrixType *Pn1        ;
    MatrixType *Hmatrix_1  ;

    MatrixType *Hmatrix     = NULL;
    MatrixType *Pn0         = NULL;
    MatrixType *Hmatrix_m1  = NULL;
    MatrixType *Hmatrix_0   = NULL;
    matrix_size = n2*sizeof(MatrixType);
    if(ct.tddft_gpu)
    {
        gpuMalloc((void **)&Hmatrix, n2*sizeof(MatrixType));
        gpuMalloc((void **)&Hmatrix_m1, n2*sizeof(MatrixType));
        gpuMalloc((void **)&Hmatrix_0, n2*sizeof(MatrixType));
        gpuMalloc((void **)&Pn0, 2*n2*sizeof(double));

        gpuMalloc((void **)&Hmatrix_1, n2*sizeof(MatrixType));
        gpuMalloc((void **)&Pn1, 2*n2*sizeof(double));
    }
    else
    {
        Pn1  = (MatrixType *)RmgMallocHost((size_t)n2*sizeof(double)*2);
    }

    if(pct.gridpe == 0) {
        printf("\n Number of states used for TDDFT: Nbasis =  %d \n",numst);
        printf(" Propagator used :  Ieldyn = %d  \\1=BCH, 2=Diagonalizer\\ \n",Ieldyn) ;
    }

    get_dipole(rho.data(), rhoc.data(), dipole_tot);
    if(ct.dipole_corr[0]+ct.dipole_corr[1]+ct.dipole_corr[2] >0)
    {
        DipoleCorrection(dipole_tot,  vh_dipole.data());
    }

    double vtxc, etxc;
    std::vector<double> Eterms_ground(6, 0.0);
    std::vector<double> Eterms_1step(6, 0.0);
    std::vector<double> Eterms(6, 0.0);
    efactor = ct.energy_output_conversion[ct.energy_output_units];
    eunits = ct.energy_output_string[ct.energy_output_units].c_str();

    RmgTimer *RT0 = new RmgTimer("2-TDDFT");

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


}

