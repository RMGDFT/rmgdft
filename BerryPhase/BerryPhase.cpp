/*
 *
 * Copyright 2023 The RMG Project Developers. See the COPYRIGHT file 
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

#include <math.h>
#include "transition.h"
#include "rmg_error.h"
#include "rmg_mangling.h"
#include "Symmetry.h"
#include "BerryPhase.h"
#include "RmgGemm.h"
#include "GlobalSums.h"
#include "blas.h"
#include "blas_driver.h"
#include "Scalapack.h"

BerryPhase *Rmg_BP;


template <typename KpointType>
void Eigen(KpointType *distA, double *eigs, KpointType *distV, int N, int M, Scalapack &Sp);

BerryPhase::BerryPhase(void)
{

    if(!ct.norm_conserving_pp)
    {
        rmg_error_handler(__FILE__, __LINE__, "only support norm-conserving pp now\n");
    }

    this->efield_mag = 0.0;
    if(ct.forceflag== TDDFT)
    {
    //    if(std::abs(ct.efield_tddft_xtal[0]) > 1.0e-10) ct.BerryPhase_dir = 0;
    //    if(std::abs(ct.efield_tddft_xtal[1]) > 1.0e-10) ct.BerryPhase_dir = 1;
    //    if(std::abs(ct.efield_tddft_xtal[2]) > 1.0e-10) ct.BerryPhase_dir = 2;
    //    this->efield_mag = ct.efield_tddft_xtal[ct.BerryPhase_dir];
    }
    else
    {
        if(std::abs(ct.efield_xtal[0]) > 1.0e-10) ct.BerryPhase_dir = 0;
        if(std::abs(ct.efield_xtal[1]) > 1.0e-10) ct.BerryPhase_dir = 1;
        if(std::abs(ct.efield_xtal[2]) > 1.0e-10) ct.BerryPhase_dir = 2;
        this->efield_mag = ct.efield_xtal[ct.BerryPhase_dir];
    }

    BerryPhase_dir = ct.BerryPhase_dir;

    int kmesh[3], is_shift[3];
    for (int i = 0; i < 3; i++)
    {
        kmesh[i] = ct.kpoint_mesh[i];
        is_shift[i] = ct.kpoint_is_shift[i];
    }
    num_kpp = ct.kpoint_mesh[BerryPhase_dir];
    if(num_kpp < 1) 
    {
        rmg_error_handler(__FILE__, __LINE__,"each string needs more than one kpoint now\n");
    }
    kmesh[BerryPhase_dir] = 1;
    is_shift[BerryPhase_dir] = 0;
    num_kort = init_kpoints(kmesh, is_shift);

    // now ct.kp only have the 2d k point orthogonal to the berry phase direction
    boost::multi_array<double, 2> kp_2d;
    kp_2d.resize(boost::extents[num_kort][3]);


    kweight_string.resize(num_kort);

    for(int i = 0; i < num_kort; i++)
    {
        kp_2d[i][0] = ct.kp[i].kpt[0];
        kp_2d[i][1] = ct.kp[i].kpt[1];
        kp_2d[i][2] = ct.kp[i].kpt[2];
        kweight_string[i] = ct.kp[i].kweight;
    }


    ct.num_kpts = num_kort * num_kpp;
    ct.kp.resize(ct.num_kpts);

    double dk = 1.0/(double)(num_kpp);
    double kpp_start = dk *0.5 * ct.kpoint_is_shift[BerryPhase_dir];
    for(int iort = 0; iort < num_kort; iort++)
    {
        for(int jpp = 0; jpp < num_kpp; jpp++)
        {
            ct.kp[iort * num_kpp + jpp].kpt[0] = kp_2d[iort][0];
            ct.kp[iort * num_kpp + jpp].kpt[1] = kp_2d[iort][1];
            ct.kp[iort * num_kpp + jpp].kpt[2] = kp_2d[iort][2];
            ct.kp[iort * num_kpp + jpp].kweight = kweight_string[iort]/(double)num_kpp;
            ct.kp[iort * num_kpp + jpp].kpt[BerryPhase_dir] = kpp_start + (jpp-num_kpp/2) * dk;
        }
    }

    for (int kpt = 0; kpt < ct.num_kpts; kpt++) {
        double v1, v2, v3;

        v1 = ct.kp[kpt].kpt[0] *Rmg_L.b0[0]
            + ct.kp[kpt].kpt[1] *Rmg_L.b1[0] 
            + ct.kp[kpt].kpt[2] *Rmg_L.b2[0];
        v2 = ct.kp[kpt].kpt[0] *Rmg_L.b0[1]
            + ct.kp[kpt].kpt[1] *Rmg_L.b1[1] 
            + ct.kp[kpt].kpt[2] *Rmg_L.b2[1];
        v3 = ct.kp[kpt].kpt[0] *Rmg_L.b0[2]
            + ct.kp[kpt].kpt[1] *Rmg_L.b1[2] 
            + ct.kp[kpt].kpt[2] *Rmg_L.b2[2];

        ct.kp[kpt].kvec[0] = v1 * twoPI;
        ct.kp[kpt].kvec[1] = v2 * twoPI;
        ct.kp[kpt].kvec[2] = v3 * twoPI;
        ct.kp[kpt].kmag = (v1 * v1 + v2 * v2 + v3 * v3) * twoPI * twoPI;
    }

    if (ct.verbose)
    {
        rmg_printf("\n num_k %d num_string(pe) %d(%d) k_in_string %d in unit of reciprocal lattice", ct.num_kpts, num_kort, num_kort_pe, num_kpp);
        for(int kpt = 0; kpt < ct.num_kpts; kpt++)
            rmg_printf("\n kvec %d  %f %f %f %f \n", kpt, ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ct.kp[kpt].kweight);
    }


}

void BerryPhase::init(Kpoint<double> **Kptr)
{
    std::cout << "not programed yet for gamma point" << std::endl;
    rmg_error_handler(__FILE__, __LINE__, "only support complex version-non-gamma now\n");
}
void BerryPhase::init(Kpoint<std::complex<double>> **Kptr)
{
    // find the number of states with non-zero occupation
    nband_occ = ct.nel/2;
    pbasis = Rmg_G->get_P0_BASIS(1);
    pbasis_noncoll = pbasis * ct.noncoll_factor;
    vel = Rmg_L.get_omega() / ((double)((size_t)Rmg_G->get_NX_GRID(1) * (size_t)Rmg_G->get_NY_GRID(1) * (size_t)Rmg_G->get_NZ_GRID(1)));
    vel_C = vel;
    wfc_size = nband_occ * pbasis_noncoll * sizeof(std::complex<double>);
    mat = (std::complex<double> *)RmgMallocHost((size_t)nband_occ * ct.max_states*sizeof(std::complex<double>));
    psi_k0 = (std::complex<double> *)RmgMallocHost(wfc_size);

    //eq 15 of PRB 47, 1651(1993) King-Smith and Vanderbilt
    // phik in eq 15
    // mat = <psi_k | psi_(k+1)>
    // BP_matrix_cpu: <psi_k |psi_(k+1)> ^-1

    if(std::abs(efield_mag) > eps )
    {
        for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {
            Kptr[kpt]->BP_matrix_cpu = (std::complex<double> *)RmgMallocHost((size_t)nband_occ * nband_occ*sizeof(std::complex<double>));
            Kptr[kpt]->BP_Gnk = (std::complex<double> *)RmgMallocHost((size_t)nband_occ * pbasis_noncoll*sizeof(std::complex<double>));
            Kptr[kpt]->BP_psi = (std::complex<double> *)RmgMallocHost((size_t)nband_occ * pbasis_noncoll*sizeof(std::complex<double>));
            for(int idx = 0; idx < nband_occ * nband_occ; idx++) Kptr[kpt]->BP_matrix_cpu [idx] = 0.0;
            for(int idx = 0; idx < nband_occ * pbasis_noncoll;  idx++) 
            {
                Kptr[kpt]->BP_Gnk [idx] = 0.0;
                Kptr[kpt]->BP_psi [idx] = 0.0;
            }
        }
    } 
}

void BerryPhase::CalcBP (Kpoint<double> **Kptr)
{
    rmg_error_handler(__FILE__, __LINE__," only support non-gamma point now\n");
}
void BerryPhase::CalcBP (Kpoint<std::complex<double>> **Kptr)
{
    // following the procedure in QE bp_c_phase.f90 


    std::complex<double> beta(0.0);

    std::complex<double> det, zeta;
    std::vector<double> phik(num_kort);
    std::vector<double> pdl_elec(num_kort);
    double pdl_elec_tot;
    std::vector<std::complex<double>> cphik(num_kort);
    int info, *ipiv;
    ipiv = new int[nband_occ];
    std::complex<double> *psi_k, *psi_k1;

    double gr[3]{0.0,0.0,0.0};
    gr[BerryPhase_dir] = -1.0;
    std::complex<double> phase;


    std::fill(phik.begin(), phik.end(), 0.0);
    std::fill(cphik.begin(), cphik.end(), 0.0);
    for(int iort = 0; iort < num_kort_pe; iort++)
    {
        int iort_gl = kort_start + iort;
        // for each string, calculate the Berry Phase
        zeta = 1.0;
        //psi_k0: the first k point in a string
        // psi_k0 * exp(i b r) for the k+1 of the last k point in a string
        // b: reciprocal lattice vector along Berry phase direction
        // if BerryPhase direction along b0, only real space along a0 direction needed
        // exp(i ix/Nx * 2PI) or exp( i iy/Ny * 2Pi), ..
        memcpy(psi_k0, Kptr[iort*num_kpp]->orbital_storage, wfc_size);

        RmgTimer *RT1 = new RmgTimer("Berry Phase: psi0 * exp(-igr)");
        psi_x_phase(psi_k0, gr, nband_occ);
        delete RT1;
        RT1 = new RmgTimer("Berry Phase: phase calculation");
        for(int jpp = 0; jpp < num_kpp; jpp++)
        {
            int ik_index = iort * num_kpp + jpp;
            psi_k = Kptr[ik_index]->orbital_storage;
            if(std::abs(efield_mag) > eps ) memcpy(Kptr[ik_index]->BP_psi, psi_k, wfc_size);
            std::complex<double> *BP_matrix_cpu = Kptr[ik_index]->BP_matrix_cpu;
            if(jpp == num_kpp-1)
            {
                psi_k1 = psi_k0;
            }
            else
            {
                psi_k1 = Kptr[ik_index+1]->orbital_storage;
            }

            RmgGemm("c", "n", nband_occ, nband_occ, pbasis_noncoll, vel_C, psi_k, pbasis_noncoll, psi_k1, pbasis_noncoll, beta, mat, nband_occ);
            BlockAllreduce(mat, nband_occ * nband_occ, pct.grid_comm);

            //calculate determinant of mat <psi_k |psi_(k+1)>

            zgetrf(&nband_occ, &nband_occ, (double *)mat, &nband_occ, ipiv, &info);
            if (info != 0)
            {
                rmg_printf ("error in zgetrf BerryPhase.cpp with INFO = %d \n", info);
                fflush (NULL);
                rmg_error_handler(__FILE__, __LINE__,"zgetrf failed\n");
            }

            det = 1.0;
            for (int i = 0; i < nband_occ; i++)
            {
                det = det * mat[i * nband_occ +i];
                if(i != ipiv[i]) det = -det;
            }


            if(std::abs(efield_mag) > eps )
            {
                int lwork = nband_occ * nband_occ;
                zgetri(&nband_occ, (double *)mat, &nband_occ, ipiv, (double *)BP_matrix_cpu, &lwork, &info);
                if (info != 0)
                {
                    rmg_printf ("error in zgetri BerryPhase.cpp with INFO = %d \n", info);
                    fflush (NULL);
                    rmg_error_handler(__FILE__, __LINE__,"zgetri failed\n");
                }
                memcpy(BP_matrix_cpu, mat, nband_occ * nband_occ * sizeof(std::complex<double>));
            }


            if(ct.verbose)
            {
                rmg_printf("kort %d kpp %d  det %f %f\n", iort, jpp, std::real(det), std::imag(det));
            }
            zeta = zeta * det;

        }
        if(std::abs(zeta) < eps)
        {

            phik[iort_gl] = 0.0;
            cphik[iort_gl] = 1.0;
        }
        else
        {
            phik[iort_gl] = std::imag( log(zeta) );
            cphik[iort_gl] = std::complex<double>(cos(phik[iort_gl]), sin(phik[iort_gl]));
            if(ct.verbose)
            {
                rmg_printf("kort %d phik  %f %f %f\n", iort,  phik[iort_gl], zeta);
            }
        }
        delete RT1;
    }

    MPI_Allreduce(MPI_IN_PLACE, phik.data(), num_kort, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
    MPI_Allreduce(MPI_IN_PLACE, cphik.data(), num_kort, MPI_DOUBLE_COMPLEX, MPI_SUM, pct.kpsub_comm);

    //  -------------------------------------------------------------------------   !
    //                    electronic polarization: phase average                    !
    //  -------------------------------------------------------------------------   !

    //  --- Initialize average of phases as complex numbers ---
    std::complex<double> cave(0.0);
    double  phik_ave(0.0);

    for(int iort = 0; iort < num_kort; iort++)
    {
        cave += kweight_string[iort] * cphik[iort];
    }

    //     --- Get the angle corresponding to the complex numbers average ---
    double theta0=atan2(std::imag(cave), std::real(cave));
    //     --- Put the phases in an around theta0 ---
    for(int iort = 0; iort < num_kort; iort++)
    {
        cphik[iort] = cphik[iort]/cave;
        double dtheta=atan2(std::imag(cphik[iort]), std::real(cphik[iort]));
        phik[iort]=theta0+dtheta;
        //rmg_printf("kort %d phik after ave %f %f\n", iort,  phik[iort]);
        // take mod so phase is -Pi to Pi
        phik[iort] = phik[iort] - PI * std::round(phik[iort]/PI); 
    }

    // you need to fix jumps before you take average
    double t1=phik[0]/PI;
    for(int iort = 0; iort < num_kort; iort++)
    {
        double t = phik[iort]/PI;
        if(abs(t+1.0-t1) < abs(t-t1))phik[iort]=phik[iort]+PI;
        if(abs(t-1.0-t1) < abs(t-t1))phik[iort]=phik[iort]-PI;
        pdl_elec[iort] = phik[iort]/PI ;
        phik_ave=phik_ave+kweight_string[iort]*phik[iort];
        if(ct.verbose)
        {
            rmg_printf("\n kstring %d weight %f phase %f", iort, kweight_string[iort], pdl_elec[iort]);
        }
    }

    pdl_elec_tot = phik_ave/twoPI;
    if(ct.nspin == 1) pdl_elec_tot *=2.0;
    // spin sum 
    pdl_elec_tot = pdl_elec_tot - 2.0*std::round(pdl_elec_tot/2.0);
    // pdl_elec_tot is [-1.0, 1.0] 

    // ionic phase  z * r_ion * b_berryPhaseDir

    double pdl_ion_tot = 0.0;
    for(int ion = 0; ion < ct.num_ions; ion++)
    {
        ION *iptr = &Atoms[ion];
        double Zj = Species[iptr->species].zvalence;
        pdl_ion_tot += Zj * iptr->xtal[BerryPhase_dir];
    }
    pdl_ion_tot = pdl_ion_tot - 2.0 * std::round(pdl_ion_tot/2.0);
    double pdl_tot = pdl_elec_tot + pdl_ion_tot;
    pdl_tot = pdl_tot - 2.0 * std::round(pdl_tot/2.0);

    rmg_printf("\n  Electronic phase %f ", pdl_elec_tot);
    rmg_printf("\n  Ionic      phase %f ", pdl_ion_tot);
    rmg_printf("\n  Total      phase %f ", pdl_tot);

    //  Polarization 

    // adapted from QE 
    //    Calculate direction of polarization and modulus of lattice vector ---
    //    lattice vector or reciprocal lattice vector with direction
    double rmod(1.0);
    if(BerryPhase_dir == 0)
    {
        rmod = Rmg_L.get_xside();
    }
    if(BerryPhase_dir == 1)
    {
        rmod = Rmg_L.get_yside();
    }
    if(BerryPhase_dir == 2)
    {
        rmod = Rmg_L.get_zside();
    }
    this->eai = rmod * this->efield_mag;
    //  --- Give polarization in units of (e/Omega).bohr ---
    pol_elec = pdl_elec_tot * rmod;
    pol_ion = pdl_ion_tot * rmod;  
    pol_tot = pdl_tot * rmod;
    rmg_printf("\n  Polarization at direction %d  = %e (e/Omega)*bohr", BerryPhase_dir, pdl_tot * rmod);
    rmg_printf("\n  Polarization at direction %d  = %e e/bohr^2", BerryPhase_dir, pdl_tot * rmod/Rmg_L.omega);
    rmg_printf("\n  Polarization at direction %d  = %e C/m^2", BerryPhase_dir, pdl_tot * rmod/Rmg_L.omega * e_C /(a0_SI * a0_SI));

    delete [] ipiv;

    if(std::abs(efield_mag) > eps )
    {
        Calc_Gnk(Kptr);
    }

}

// calculate second part of eq. 4  
//  Ivo Souza, Jorge I´n˜iguez, and David Vanderbilt, PRL2002, 117602
// S^-1(k, k+1) psi_(k+1) - S^1(k, k-1) psi_(k-1)
// S^-1(k, k-1) = [S^-1(k-1,k)]^+ hermitian
//

void BerryPhase::Calc_Gnk (Kpoint<std::complex<double>> **Kptr)
{
    std::complex<double> alpha = std::complex<double>(0.0, -this->eai/twoPI/2.0 * num_kpp);
    //rmg_printf("\n eai alpha %f %f", alpha);
    std::complex<double> zero(0.0), mone(-1.0);
    std::complex<double> *psi_kp1=NULL, *psi_km1=NULL;
    std::complex<double> *BP_matrix_p1 =NULL, *BP_matrix_m1;
    std::complex<double> *Gnk =NULL;
    double gr[3]{0.0, 0.0, 0.0};
    for(int iort = 0; iort < num_kort_pe; iort++)
    {
        for(int jpp = 0; jpp < num_kpp; jpp++)
        {
            int ik_index = iort * num_kpp + jpp;
            BP_matrix_p1 = Kptr[ik_index]->BP_matrix_cpu;
            Gnk = Kptr[ik_index]->BP_Gnk;

            BP_matrix_p1 = Kptr[iort * num_kpp + jpp]->BP_matrix_cpu;
            if(jpp == 0 )
            {
                BP_matrix_m1 = Kptr[iort * num_kpp + num_kpp-1]->BP_matrix_cpu;
                memcpy(psi_k0, Kptr[iort*num_kpp + num_kpp-1]->orbital_storage, wfc_size);
                gr[BerryPhase_dir] = 1.0;
                psi_x_phase(psi_k0, gr, nband_occ);
                psi_km1 = psi_k0;
                psi_kp1 = Kptr[iort*num_kpp + + jpp +1 ]->orbital_storage;
            }
            else if(jpp == num_kpp -1)
            {
                BP_matrix_m1 = Kptr[iort * num_kpp + jpp-1]->BP_matrix_cpu;
                memcpy(psi_k0, Kptr[iort*num_kpp ]->orbital_storage, wfc_size);
                gr[BerryPhase_dir] = -1.0;
                psi_x_phase(psi_k0, gr, nband_occ);
                psi_km1 = Kptr[iort*num_kpp + + jpp -1 ]->orbital_storage;
                psi_kp1 = psi_k0;
            }
            else
            {
                BP_matrix_m1 = Kptr[iort * num_kpp + jpp-1]->BP_matrix_cpu;
                psi_km1 = Kptr[iort*num_kpp + + jpp -1 ]->orbital_storage;
                psi_kp1 = Kptr[iort*num_kpp + + jpp +1 ]->orbital_storage;
            }

            RmgGemm("N", "N", pbasis_noncoll, nband_occ, nband_occ, alpha,
                    psi_kp1, pbasis_noncoll, BP_matrix_p1, nband_occ, 
                    zero, Gnk, pbasis_noncoll);
            RmgGemm("N", "C", pbasis_noncoll, nband_occ, nband_occ, alpha,
                    psi_km1, pbasis_noncoll, BP_matrix_m1, nband_occ, 
                    mone, Gnk, pbasis_noncoll);

        }
    }

}

void BerryPhase::Apply_BP_Hpsi(Kpoint<double> *kptr, int num_states, double *psi, double *h_psi)
{
}
void BerryPhase::Apply_BP_Hpsi(Kpoint<std::complex<double>> *kptr, int num_states, std::complex<double> *psi, std::complex<double> *h_psi)
{
    if(std::abs(efield_mag) < eps ) return;

    // extra Hamiltonian operator |Gnk ><psi_bp| + |psi_bp >< Gnk|, applied to all states including unoccupied 
    // BP_Hpsi =     |Gnk > <psi_bp | psi current> 
    std::complex<double> one(1.0), zero(0.0); 
    RmgGemm("c", "n", nband_occ, num_states, pbasis_noncoll, vel_C, kptr->BP_psi, pbasis_noncoll, psi, pbasis_noncoll, zero, mat, nband_occ);
    BlockAllreduce(mat, (size_t)nband_occ * num_states, pct.grid_comm);
    RmgGemm("N", "N", pbasis_noncoll, num_states, nband_occ, one,
            kptr->BP_Gnk, pbasis_noncoll, mat, nband_occ, 
            one, h_psi, pbasis_noncoll);
    // BP_Hpsi +=     |psi_bp > <Gnk | psi current> 
    RmgGemm("c", "n", nband_occ, num_states, pbasis_noncoll, vel_C, kptr->BP_Gnk, pbasis_noncoll, psi, pbasis_noncoll, zero, mat, nband_occ);
    BlockAllreduce(mat, (size_t)nband_occ * num_states, pct.grid_comm);
    RmgGemm("N", "N", pbasis_noncoll, num_states, nband_occ, one,
            kptr->BP_psi, pbasis_noncoll, mat, nband_occ, 
            one, h_psi, pbasis_noncoll);
}

BerryPhase::~BerryPhase(void)
{
}
void BerryPhase::psi_x_phase(std::complex<double> *psi_k0, double gr[3], int nband)
{

    // psi_k0 = psi_k0 * exp(i gr )
    int pbasis = Rmg_G->get_P0_BASIS(1);
    int pbasis_noncoll = pbasis * ct.noncoll_factor;
    int px0_grid = Rmg_G->get_PX0_GRID(1);
    int py0_grid = Rmg_G->get_PY0_GRID(1);
    int pz0_grid = Rmg_G->get_PZ0_GRID(1);
    int pxoffset = Rmg_G->get_PX_OFFSET(1);
    int pyoffset = Rmg_G->get_PY_OFFSET(1);
    int pzoffset = Rmg_G->get_PZ_OFFSET(1);
    int Nx = Rmg_G->get_NX_GRID(1);
    int Ny = Rmg_G->get_NY_GRID(1);
    int Nz = Rmg_G->get_NZ_GRID(1);

    for(int st = 0; st < nband; st++)
    {
        for (int ix = 0; ix < px0_grid; ix++)
        {
            for(int iy = 0; iy < py0_grid; iy++)
            {
                for(int iz = 0; iz < pz0_grid; iz++)
                {
                    int idx = ix * py0_grid * pz0_grid + iy * pz0_grid + iz;
                    double gmr = ( (ix + pxoffset)/(double)Nx * gr[0] +
                            (iy + pyoffset)/(double)Ny * gr[1] +
                            (iz + pzoffset)/(double)Nz * gr[2] ) * twoPI;
                    psi_k0[st * pbasis_noncoll + idx] *= std::complex(cos(gmr), sin(gmr));
                }
            }
        }

    }
}

void BerryPhase::CalcBP_Skk1 (Kpoint<double> **Kptr, int tddft_start_state, double *matrix_glob, Scalapack &Sp )
{
    rmg_error_handler(__FILE__, __LINE__," only support non-gamma point now\n");
}
void BerryPhase::CalcBP_Skk1 (Kpoint<std::complex<double>> **Kptr, int tddft_start_state, std::complex<double> *mat_glob, Scalapack &Sp )
{

    //calculating <Psi_k | Psi_k+1>
    std::complex<double> beta(0.0);

    std::complex<double> *psi_k, *psi_k1;

    double gr[3]{0.0,0.0,0.0};
    gr[BerryPhase_dir] = -1.0;
    std::complex<double> phase;

    int numst = ct.num_states - tddft_start_state;
    if(numst != Sp.GetN())
    {
        rmg_printf("\n Scalpapack wrong !! numst = %d  N = %d \n", numst, Sp.GetN() );
        rmg_error_handler(__FILE__, __LINE__," scalapack wrong\n");
    }

    wfc_size = numst * pbasis_noncoll * sizeof(std::complex<double>);
    if(psi_k0 != NULL) RmgFreeHost(psi_k0);
    psi_k0 = (std::complex<double> *)RmgMallocHost(wfc_size);

    int Mdim = Sp.GetDistMdim();
    int Ndim = Sp.GetDistNdim();
    int n2 = Mdim * Ndim;
    int *desca = Sp.GetDistDesca();

    if(Kptr[0]->BP_Skk1_cpu == NULL)
    {
        for(int ik = 0; ik < ct.num_kpts_pe; ik++)
        {
            Kptr[ik]->BP_Skk1_cpu
                 = (std::complex<double> *)RmgMallocHost((size_t)n2*sizeof(std::complex<double>));

        }
    }

    
    RmgTimer *RT1 = new RmgTimer("Berry Phase: tddft Skk1");
    for(int iort = 0; iort < num_kort_pe; iort++)
    {
        //psi_k0: the first k point in a string
        // psi_k0 * exp(i b r) for the k+1 of the last k point in a string
        // b: reciprocal lattice vector along Berry phase direction
        // if BerryPhase direction along b0, only real space along a0 direction needed
        // exp(i ix/Nx * 2PI) or exp( i iy/Ny * 2Pi), ..
        memcpy(psi_k0, Kptr[iort*num_kpp]->orbital_storage + tddft_start_state * pbasis_noncoll, wfc_size);

        psi_x_phase(psi_k0, gr, numst);
        for(int jpp = 0; jpp < num_kpp; jpp++)
        {
            int ik_index = iort * num_kpp + jpp;
            psi_k = Kptr[ik_index]->orbital_storage + tddft_start_state * pbasis_noncoll;
            if(jpp == num_kpp-1)
            {
                psi_k1 = psi_k0;
            }
            else
            {
                psi_k1 = Kptr[ik_index+1]->orbital_storage;
            }

            RmgGemm("c", "n", numst, numst, pbasis_noncoll, vel_C, psi_k, pbasis_noncoll, psi_k1, pbasis_noncoll, beta, mat_glob, numst);
            BlockAllreduce(mat_glob, numst * numst, pct.grid_comm);

            Sp.CopySquareMatrixToDistArray(mat_glob, Kptr[ik_index]->BP_Skk1_cpu, numst, desca);

        }
    }
    delete RT1;

}
void BerryPhase::CalcBP_tddft (Kpoint<double> **Kptr, double &tot_bp_pol, double *mat_glob, Scalapack &Sp)
{
    rmg_error_handler(__FILE__, __LINE__," only support non-gamma point now\n");
}
void BerryPhase::CalcBP_tddft (Kpoint<std::complex<double>> **Kptr, double &tot_bp_pol, std::complex<double> *mat_glob, Scalapack &Sp)
{



    std::complex<double> det, zeta;
    std::vector<double> phik(num_kort);
    std::vector<double> pdl_elec(num_kort);
    double pdl_elec_tot;
    std::vector<std::complex<double>> cphik(num_kort);
    int info, *ipiv;
    ipiv = new int[nband_occ];

    int nprow = Sp.GetRows();
    int npcol = Sp.GetCols();
    int NB = Sp.GetNB();
    int *desca = Sp.GetDistDesca();
    int numst = Sp.GetN();
    int ione = 1;

    int Mdim = Sp.GetDistMdim();
    int Ndim = Sp.GetDistNdim();
    int n2 = Sp.GetDistMdim() * Sp.GetDistNdim();


    double *eigs = new double[numst];
    std::complex<double> *Cmat = new std::complex<double>[n2];
    std::complex<double> *CijSkk1 = new std::complex<double>[n2];
    for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        // Pn1_cpu will store the eigenvectors
        // new psi = Pn1_cpu * psi
        Eigen(Kptr[kpt]->Hmatrix_1_cpu, eigs, Kptr[kpt]->Pn1_cpu, numst, numst, Sp);  

    }

    std::complex<double> alpha(1.0), beta(0.0);
    RmgTimer  *RT1 = new RmgTimer("Berry Phase: tddft phase calculation");
    std::fill(phik.begin(), phik.end(), 0.0);
    std::fill(cphik.begin(), cphik.end(), 0.0);
    for(int iort = 0; iort < num_kort_pe; iort++)
    {
        int iort_gl = kort_start + iort;
        // for each string, calculate the Berry Phase
        // BP_Skk1 = <psi_k |psi_k1>
        zeta = 1.0;

        for(int jpp = 0; jpp < num_kpp; jpp++)
        {
            int ik_index = iort * num_kpp + jpp;

            std::complex<double> *Cij_k = Kptr[ik_index]->Pn1_cpu;
            std::complex<double> *Cij_k1 = NULL;
            if(jpp == num_kpp-1)
            {
                Cij_k1= Kptr[iort * num_kpp]->Pn1_cpu;
            }
            else
            {
                Cij_k1= Kptr[ik_index+1]->Pn1_cpu;
            } 
            std::complex<double> *Skk1 = Kptr[ik_index]->BP_Skk1_cpu;
            zgemm_driver ("C", "N", numst, numst, numst, alpha, Cij_k, ione, ione, desca,
                    Skk1, ione, ione, desca, beta, CijSkk1, ione, ione, desca);
            zgemm_driver ("N", "N", numst, numst, numst, alpha, CijSkk1, ione, ione, desca,
                    Cij_k1, ione, ione, desca, beta, Cmat, ione, ione, desca);

            Sp.GatherEigvectors(mat_glob, Cmat);

            zgetrf(&nband_occ, &nband_occ, (double *)mat_glob, &numst, ipiv, &info);
            if (info != 0)
            {
                rmg_printf ("error in zgetrf BerryPhase.cpp with INFO = %d \n", info);
                fflush (NULL);
                rmg_error_handler(__FILE__, __LINE__,"zgetrf failed\n");
            }

            det = 1.0;
            for (int i = 0; i < nband_occ; i++)
            {
                det = det * mat_glob[i * numst +i];
                if(i != ipiv[i]) det = -det;
            }


            if(ct.verbose)
            {
                rmg_printf("kort %d kpp %d  det %f %f\n", iort, jpp, std::real(det), std::imag(det));
            }
            zeta = zeta * det;

        }
        if(std::abs(zeta) < eps)
        {

            phik[iort_gl] = 0.0;
            cphik[iort_gl] = 1.0;
        }
        else
        {
            phik[iort_gl] = std::imag( log(zeta) );
            cphik[iort_gl] = std::complex<double>(cos(phik[iort_gl]), sin(phik[iort_gl]));
            if(ct.verbose)
            {
                rmg_printf("kort %d phik  %f %f %f\n", iort,  phik[iort_gl], zeta);
            }
        }
    }
    delete RT1;

    MPI_Allreduce(MPI_IN_PLACE, phik.data(), num_kort, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
    MPI_Allreduce(MPI_IN_PLACE, cphik.data(), num_kort, MPI_DOUBLE_COMPLEX, MPI_SUM, pct.kpsub_comm);

    //  -------------------------------------------------------------------------   !
    //                    electronic polarization: phase average                    !
    //  -------------------------------------------------------------------------   !

    //  --- Initialize average of phases as complex numbers ---
    std::complex<double> cave(0.0);
    double  phik_ave(0.0);

    for(int iort = 0; iort < num_kort; iort++)
    {
        cave += kweight_string[iort] * cphik[iort];
    }

    //     --- Get the angle corresponding to the complex numbers average ---
    double theta0=atan2(std::imag(cave), std::real(cave));
    //     --- Put the phases in an around theta0 ---
    for(int iort = 0; iort < num_kort; iort++)
    {
        cphik[iort] = cphik[iort]/cave;
        double dtheta=atan2(std::imag(cphik[iort]), std::real(cphik[iort]));
        phik[iort]=theta0+dtheta;
        //rmg_printf("kort %d phik after ave %f %f\n", iort,  phik[iort]);
        // take mod so phase is -Pi to Pi
        phik[iort] = phik[iort] - PI * std::round(phik[iort]/PI); 
    }

    // you need to fix jumps before you take average
    double t1=phik[0]/PI;
    for(int iort = 0; iort < num_kort; iort++)
    {
        double t = phik[iort]/PI;
        if(abs(t+1.0-t1) < abs(t-t1))phik[iort]=phik[iort]+PI;
        if(abs(t-1.0-t1) < abs(t-t1))phik[iort]=phik[iort]-PI;
        pdl_elec[iort] = phik[iort]/PI ;
        phik_ave=phik_ave+kweight_string[iort]*phik[iort];
        if(ct.verbose)
        {
            rmg_printf("\n kstring %d weight %f phase %f", iort, kweight_string[iort], pdl_elec[iort]);
        }
    }

    pdl_elec_tot = phik_ave/twoPI;
    if(ct.nspin == 1) pdl_elec_tot *=2.0;
    // spin sum 
    pdl_elec_tot = pdl_elec_tot - 2.0*std::round(pdl_elec_tot/2.0);
    // pdl_elec_tot is [-1.0, 1.0] 

    // ionic phase  z * r_ion * b_berryPhaseDir

    double pdl_ion_tot = 0.0;
    for(int ion = 0; ion < ct.num_ions; ion++)
    {
        ION *iptr = &Atoms[ion];
        double Zj = Species[iptr->species].zvalence;
        pdl_ion_tot += Zj * iptr->xtal[BerryPhase_dir];
    }
    pdl_ion_tot = pdl_ion_tot - 2.0 * std::round(pdl_ion_tot/2.0);
    double pdl_tot = pdl_elec_tot + pdl_ion_tot;
    pdl_tot = pdl_tot - 2.0 * std::round(pdl_tot/2.0);

    if(ct.verbose)
    {
        rmg_printf("\n  Electronic phase %f ", pdl_elec_tot);
        rmg_printf("\n  Ionic      phase %f ", pdl_ion_tot);
        rmg_printf("\n  Total      phase %f ", pdl_tot);
    }

    //  Polarization 

    // adapted from QE 
    //    Calculate direction of polarization and modulus of lattice vector ---
    //    lattice vector or reciprocal lattice vector with direction
    double rmod(1.0);
    if(BerryPhase_dir == 0)
    {
        rmod = Rmg_L.get_xside();
    }
    if(BerryPhase_dir == 1)
    {
        rmod = Rmg_L.get_yside();
    }
    if(BerryPhase_dir == 2)
    {
        rmod = Rmg_L.get_zside();
    }
    this->eai = rmod * this->efield_mag;
    //  --- Give polarization in units of (e/Omega).bohr ---
    pol_elec = pdl_elec_tot * rmod;
    pol_ion = pdl_ion_tot * rmod;  
    pol_tot = pdl_tot * rmod;
    if(ct.verbose)
    {
        rmg_printf("\n  Polarization at direction %d  = %e (e/Omega)*bohr", BerryPhase_dir, pdl_tot * rmod);
        rmg_printf("\n  Polarization at direction %d  = %e e/bohr^2", BerryPhase_dir, pdl_tot * rmod/Rmg_L.omega);
        rmg_printf("\n  Polarization at direction %d  = %e C/m^2", BerryPhase_dir, pdl_tot * rmod/Rmg_L.omega * e_C /(a0_SI * a0_SI));
    }

    tot_bp_pol = pdl_tot * rmod/Rmg_L.omega * e_C /(a0_SI * a0_SI);
    delete [] ipiv;
    delete [] eigs;
    delete [] Cmat;
    delete [] CijSkk1;


}

