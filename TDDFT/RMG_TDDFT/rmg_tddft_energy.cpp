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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <filesystem>
#if !(defined(_WIN32) || defined(_WIN64))
#include <unistd.h>
#else
#include <io.h>
#endif


#include <limits.h>
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "rmg_reduce.h"
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
#include "rmg_sum_all.h"
#include "GatherScatter.h"
#include "rmg_hvector.h"
#include "blacs.h"
#include "rmg_tddft.h"



template void rmg::tddft<double, double>::TddftEnergyInit ( spinobj<double> &vxc,
        fgobj<double> &vh, fgobj<double> &vnuc,
        spinobj<double> &rho_ground,
        fgobj<double> &rhocore, fgobj<double> &rhoc,
        Kpoint<double> **Kptr, 
        Scalapack &SP, int Mdim, int Ndim, std::vector<double> &Eterms);

template void rmg::tddft<double, std::complex<double>>::TddftEnergyInit ( spinobj<double> &vxc,
        fgobj<double> &vh, fgobj<double> &vnuc,
        spinobj<double> &rho_ground,
        fgobj<double> &rhocore, fgobj<double> &rhoc, 
        Kpoint<double> **, 
        Scalapack &SP, int Mdim, int Ndim, std::vector<double> &Eterms);

template void rmg::tddft<std::complex<double>, std::complex<double>>::TddftEnergyInit ( spinobj<double> &vxc,
        fgobj<double> &vh, fgobj<double> &vnuc,
        spinobj<double> &rho_ground,
        fgobj<double> &rhocore, fgobj<double> &rhoc, 
        Kpoint<std::complex<double>> **, 
        Scalapack &SP, int Mdim, int Ndim, std::vector<double> &Eterms);



template <typename OrbitalType, typename MatrixType> void rmg::tddft<OrbitalType, MatrixType>::TddftEnergyInit ( spinobj<double> &vxc,
        fgobj<double> &vh, fgobj<double> &vnuc,
        spinobj<double> &rho_ground,
        fgobj<double> &rhocore, fgobj<double> &rhoc,
        Kpoint<OrbitalType> **Kptr, 
        Scalapack &SP, int Mdim, int Ndim, std::vector<double> &Eterms)
{


    double vel = get_vel_f();
    std::string filename;
    int numst;

    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    numst = ct.num_states - ct.tddft_start_state; 

    double ES_0 = 0.0, ES_0a = 0.0, ES_0b = 0.0; 
    double EkinPseudo_0 = 0.0, totalE_0=0.0;
    double E_downfold = 0.0;
    double vtxc, etxc, etxc_0=0.0;

    // calculate the total energy terms of ground states
    // and Hcore matrix

    RmgTimer *RT1 = new RmgTimer("2-TDDFT: exchange/correlation");
    compute_vxc(rho_ground.data(), rhocore.data(), etxc, vtxc, vxc.data(), ct.nspin);
    etxc_0 = etxc;
    delete RT1;

    RT1 = new RmgTimer("2-TDDFT: Vh");
    VhDriver(rho_ground.data(), rhoc.data(), vh.data(), ct.vh_ext, 1.0-12);
    delete RT1;

    ct.II = IonIonEnergy_Ewald();
    // vh -> rho-rhoc
    // ES_0a:  (rho + rhoc ) * (rho-rhoc)
    // ES_0b:  (rho  * (rho-rhoc)
    // ct.ES_rhoc:  rhoc * rhoc/(r-r)
    // canceling part  from vnuc + compensating charge potential
    if (ct.nspin==2)
    {   
        /* Add the compensating charge to total charge to calculation electrostatic energy */
        for (int idx = 0; idx < FP0_BASIS; idx++) 
        {
            ES_0a += (rho_ground.up[idx] + rho_ground.dw[idx] + rhoc[idx]) * vh[idx];
            ES_0b += (rho_ground.up[idx] + rho_ground.dw[idx] ) * vh[idx];
        }

    }
    else
    {
        for (int idx = 0; idx < FP0_BASIS; idx++) 
        {
            ES_0a += (rho_ground[idx] + rhoc[idx]) * vh[idx];
            ES_0b += rho_ground[idx] * vh[idx];
        }
    }

    ES_0a = 0.5 * vel * rmg::sum_all(ES_0a, pct.grid_comm);
    ES_0b =  vel * rmg::sum_all(ES_0b, pct.grid_comm);

    ES_0 = ES_0b - ES_0a + ct.ES_rhoc;
    int nstates = Kptr[0]->nstates;
    rmg::hvector<OrbitalType> Hcore(nstates * nstates);

    wfobj<double> vtot_psi;
    wf_spinobj<double> vxc_psi;
    GetVtotPsi (vtot_psi.data(), vnuc.data(), Rmg_G->default_FG_RATIO);
    EkinPseudo_0 = 0.0;
    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {

        Kptr[kpt]->Hcore_tddft.resize(Mdim * Ndim);
        Kptr[kpt]->ComputeHcore(vtot_psi.data(), vxc_psi.data(), Hcore.data(), NULL, NULL);

        for (int st = 0; st < ct.tddft_start_state; st++)
        {
            double occ = Kptr[kpt]->Kstates[st].occupation[0] * Kptr[kpt]->kp.kweight;
            E_downfold += std::real(Hcore[st * nstates + st]) * occ;
        }

        for (int st = 0; st < numst; st++)
        {
            double occ = Kptr[kpt]->Kstates[st + ct.tddft_start_state].occupation[0] * Kptr[kpt]->kp.kweight;
            EkinPseudo_0 += std::real(Hcore[(st + ct.tddft_start_state) * nstates + st + ct.tddft_start_state]) * occ;
        }

        if(ct.tddft_tiledMM)
        {
            int numst_pe = numst/pct.local_comm_npes;

            for(int st1 = 0; st1 < numst_pe; st1++) {
                for(int st2 = 0; st2 < numst; st2++) {
                    int st1g = st1 + numst_pe * pct.local_rank;
                    Kptr[kpt]->Hcore_tddft[st1 * numst + st2] = Hcore[(st1g+ct.tddft_start_state) * nstates + st2 + ct.tddft_start_state];
                }
            }
        }
        else
        {
            int *desca = SP.GetDistDesca();
            int ictxt=desca[1], mb=desca[4], nb=desca[5], mxllda = desca[8];
            int mycol, myrow, nprow, npcol;
            Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
            int izero = 0;
            for(int i = 0; i < SP.GetDistMdim(); i++)
            {
                int i1 = i+1;
                int i_glob = indxl2g(&i1, &mb, &myrow, &izero, &nprow) + ct.tddft_start_state;
                i_glob -=1;
                for(int j = 0; j < SP.GetDistNdim(); j++)
                {
                    int j1 = j+1;
                    int j_glob = indxl2g(&j1, &nb, &mycol, &izero, &npcol) + ct.tddft_start_state;
                    j_glob -=1;
                    Kptr[kpt]->Hcore_tddft[i + j * mxllda] = Hcore[i_glob + j_glob * nstates];
                }
            }

        }
    }

    E_downfold = rmg::sum_all(E_downfold, pct.kpsub_comm);
    E_downfold = rmg::sum_all(E_downfold, pct.spin_comm);
    EkinPseudo_0 = rmg::sum_all(EkinPseudo_0, pct.kpsub_comm);
    EkinPseudo_0 = rmg::sum_all(EkinPseudo_0, pct.spin_comm);

    totalE_0 = E_downfold + EkinPseudo_0 + ES_0 + etxc_0 + ct.II;

    Eterms[0] = totalE_0     ;
    Eterms[1] = EkinPseudo_0 ;
    Eterms[2] = ES_0         ;
    Eterms[3] = etxc_0       ;
    Eterms[4] = ct.II        ;
    Eterms[5] = E_downfold   ;


}


template void rmg::tddft<double, double>::TddftEnergy ( fgobj<double> &vh, spinobj<double> &rho,fgobj<double> &rhoc,
        Kpoint<double> **Kptr, int Mdim, int Ndim, std::vector<double> &Eterms, MPI_Comm mat_comm);

template void rmg::tddft<double, std::complex<double>>::TddftEnergy (
         fgobj<double> &vh,
         spinobj<double> &rho,
         fgobj<double> &rhoc,
         Kpoint<double> **Kptr,
         int Mdim,
         int Ndim, std::vector<double> &Eterms,
         MPI_Comm mat_comm);

template void rmg::tddft<std::complex<double>, std::complex<double>>::TddftEnergy(
         fgobj<double> &vh,
         spinobj<double> &rho,
         fgobj<double> &rhoc,
         Kpoint<std::complex<double>> **Kptr,
         int Mdim,
         int Ndim, std::vector<double> &Eterms,
         MPI_Comm mat_comm);

template <typename OrbitalType, typename MatrixType> void rmg::tddft<OrbitalType, MatrixType>::TddftEnergy (
        fgobj<double> &vh,
        spinobj<double> &rho,
        fgobj<double> &rhoc,
        Kpoint<OrbitalType> **Kptr,
        int Mdim,
        int Ndim, std::vector<double> &Eterms,
        MPI_Comm mat_comm)
{


    double vel = get_vel_f();

    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    double ES = 0.0, EkinPseudo = 0.0, totalE=0.0;
    double ES_0a = 0.0, ES_0b = 0.0;
       // vh -> rho-rhoc
    // ES_0a:  (rho + rhoc ) * (rho-rhoc)
    // ES_0b:  (rho  * (rho-rhoc)
    // ct.ES_rhoc:  rhoc * rhoc/(r-r)
    // canceling part  from vnuc + compensating charge potential
    if (ct.nspin==2)
    {
        /* Add the compensating charge to total charge to calculation electrostatic energy */
        for (int idx = 0; idx < FP0_BASIS; idx++)
        {
            ES_0a += (rho.up[idx] + rho.dw[idx] + rhoc[idx]) * vh[idx];
            ES_0b += (rho.up[idx] + rho.dw[idx] ) * vh[idx];
        }

    }
    else
    {
        for (int idx = 0; idx < FP0_BASIS; idx++)
        {
            ES_0a += (rho[idx] + rhoc[idx]) * vh[idx];
            ES_0b += rho[idx] * vh[idx];
        }
    }

    ES_0a = 0.5 * vel * rmg::sum_all(ES_0a, pct.grid_comm);
    ES_0b =  vel * rmg::sum_all(ES_0b, pct.grid_comm);

    ES = ES_0b - ES_0a + ct.ES_rhoc;


    int n2 = Mdim * Ndim, ione = 1, itwo = 2;
    int n22 = n2 * sizeof(OrbitalType)/sizeof(double);
    if(typeid(OrbitalType) == typeid(MatrixType))
    {
        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) 
        {
            double kweight = Kptr[kpt]->kp.kweight;
            double dd;
            dd = ddot(&n22, (double *)Kptr[kpt]->Pn1_cpu, &ione, (double *)Kptr[kpt]->Hcore_tddft.data(), &ione);
            EkinPseudo += dd*kweight;
        }
    }
    else 
    {
        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) 
        {
            double kweight = Kptr[kpt]->kp.kweight;
            double dd;
            dd = ddot(&n2, (double *)Kptr[kpt]->Pn0_cpu, &itwo, (double *)Kptr[kpt]->Hcore_tddft.data(), &ione);
            EkinPseudo += dd*kweight;
        }
    }

    EkinPseudo = rmg::sum_all(EkinPseudo, mat_comm);
    EkinPseudo = rmg::sum_all(EkinPseudo, pct.kpsub_comm);
    EkinPseudo = rmg::sum_all(EkinPseudo, pct.spin_comm);

    Eterms[1] = EkinPseudo ;
    Eterms[2] = ES;
    totalE = 0.0;
    for(int i = 1; i < 6; i++) totalE += Eterms[i];
    Eterms[0] = totalE     ;


}
