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

BerryPhase *Rmg_BP;


BerryPhase::BerryPhase(void)
{

    if(pct.pe_kpoint != 1) 
    {
        rmg_error_handler(__FILE__, __LINE__,"no kpoint parallel now\n");
    }
    if(!ct.norm_conserving_pp)
    {
        rmg_error_handler(__FILE__, __LINE__, "only support norm-conserving pp now\n");
    }

    if(ct.forceflag== TDDFT)
    {
        if(std::abs(ct.efield_tddft_xtal[0]) > 1.0e-10) ct.BerryPhase_dir = 0;
        if(std::abs(ct.efield_tddft_xtal[1]) > 1.0e-10) ct.BerryPhase_dir = 1;
        if(std::abs(ct.efield_tddft_xtal[2]) > 1.0e-10) ct.BerryPhase_dir = 2;
        this->efield_mag = ct.efield_tddft_xtal[ct.BerryPhase_dir];
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
        rmg_printf("\n num_k %d num_string %d k_in_string %d in unit of reciprocal lattice", ct.num_kpts, num_kort, num_kpp);
        for(int kpt = 0; kpt < ct.num_kpts; kpt++)
            rmg_printf("\n kvec %d  %f %f %f %f \n", kpt, ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ct.kp[kpt].kweight);
    }



}

void BerryPhase::CalcBP (Kpoint<double> **Kptr)
{
    rmg_error_handler(__FILE__, __LINE__," only support non-gamma point now\n");
}
void BerryPhase::CalcBP (Kpoint<std::complex<double>> **Kptr)
{
    // following the procedure in QE bp_c_phase.f90 
    double eps = 1.0e-6;
    // find the number of states with non-zero occupation
    int nband_occ = 0;
    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        int st = 0;
        for(st = 0;  st < ct.num_states; st++)
        {
            if(std::abs(Kptr[kpt]->Kstates[st].occupation[0]) < eps)
            {
                break;
            }
        }
        nband_occ = std::max(nband_occ, st);
    }

    //eq 15 of PRB 47, 1651(1993) King-Smith and Vanderbilt
    // phik in eq 15
    // mat = <psi_k | psi_(k+1)>
    // BP_matrix_cpu: <psi_k |psi_(k+1)> ^-1
     
    if(std::abs(efield_mag) > 0.0 && Kptr[0]->BP_matrix_cpu == NULL)
    {
        for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {
            Kptr[kpt]->BP_matrix_cpu = (std::complex<double> *)RmgMallocHost((size_t)nband_occ * nband_occ*sizeof(std::complex<double>));
        }
    } 
    static std::complex<double> *mat=NULL, *psi_k0 = NULL;

    double vel = Rmg_L.get_omega() / ((double)((size_t)Rmg_G->get_NX_GRID(1) * (size_t)Rmg_G->get_NY_GRID(1) * (size_t)Rmg_G->get_NZ_GRID(1)));
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


    if(!mat)
    {
        mat = (std::complex<double> *)RmgMallocHost((size_t)nband_occ * nband_occ*sizeof(std::complex<double>));
        psi_k0 = (std::complex<double> *)RmgMallocHost((size_t)nband_occ * pbasis_noncoll *sizeof(std::complex<double>));
    }

    size_t wfc_size = nband_occ * pbasis_noncoll * sizeof(std::complex<double>);

    std::complex<double> alphavel(vel);
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
    gr[BerryPhase_dir] = 1.0;
    std::complex<double> phase;
    for(int iort = 0; iort < num_kort; iort++)
    {
        // for each string, calculate the Berry Phase
        zeta = 1.0;
        //psi_k0: the first k point in a string
        // psi_k0 * exp(i b r) for the k+1 of the last k point in a string
        // b: reciprocal lattice vector along Berry phase direction
        // if BerryPhase direction along b0, only real space along a0 direction needed
        // exp(i ix/Nx * 2PI) or exp( i iy/Ny * 2Pi), ..
        memcpy(psi_k0, Kptr[iort*num_kpp]->orbital_storage, wfc_size);

        RmgTimer *RT1 = new RmgTimer("Berry Phase: psi0 * exp(-igr)");
        for(int st = 0; st < nband_occ; st++)
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
                        psi_k0[st * pbasis_noncoll + idx] *= std::complex(cos(gmr), -sin(gmr));
                    }
                }
            }

        }
        delete RT1;
        for(int jpp = 0; jpp < num_kpp; jpp++)
        {
            int ik_index = iort * num_kpp + jpp;
            psi_k = Kptr[ik_index]->orbital_storage;
            std::complex<double> *BP_matrix_cpu = Kptr[ik_index]->BP_matrix_cpu;
            if(jpp == num_kpp-1)
            {
                psi_k1 = psi_k0;
            }
            else
            {
                psi_k1 = Kptr[ik_index+1]->orbital_storage;
            }

            RmgGemm("c", "n", nband_occ, nband_occ, pbasis_noncoll, alphavel, psi_k, pbasis_noncoll, psi_k1, pbasis_noncoll, beta, mat, nband_occ);
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


            if(std::abs(efield_mag) > 0.0 )
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

            phik[iort] = 0.0;
            cphik[iort] = 1.0;
        }
        else
        {
            phik[iort] = std::imag( log(zeta) );
            cphik[iort] = std::complex<double>(cos(phik[iort]), sin(phik[iort]));
            if(ct.verbose)
            {
                rmg_printf("kort %d phik  %f \n", iort,  phik[iort]);
            }
        }
    }

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
    //  --- Give polarization in units of (e/Omega).bohr ---
    pol_elec = pdl_elec_tot * rmod;
    pol_ion = pdl_ion_tot * rmod;  
    pol_tot = pdl_tot * rmod;
    rmg_printf("\n  Polarization at direction %d  = %e (e/Omega)*bohr", BerryPhase_dir, pdl_tot * rmod);
    rmg_printf("\n  Polarization at direction %d  = %e e/bohr^2", BerryPhase_dir, pdl_tot * rmod/Rmg_L.omega);
    rmg_printf("\n  Polarization at direction %d  = %e C/m^2", BerryPhase_dir, pdl_tot * rmod/Rmg_L.omega * e_C /(a0_SI * a0_SI));

    RmgFreeHost(mat);
    delete [] ipiv;
}


BerryPhase::~BerryPhase(void)
{
}

