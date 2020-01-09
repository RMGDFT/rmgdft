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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "transition.h"
#include "Atomic.h"
#include "RmgParallelFft.h"
#include "prototypes_rmg.h"
#include "GlobalSums.h"
#include <boost/math/special_functions/erf.hpp>
#include "GpuAlloc.h"
#include "RmgGemm.h"
#include "transition.h"
#include "Stress.h"
#include "AtomicInterpolate.h"
#include "RmgException.h"
#include "Functional.h"
#include "blas.h"

#include "rmg_mangling.h"
#define get_inlc                RMG_FC_MODULE(funct,get_inlc,mod_FUNCT,GET_INLC)

extern "C" int get_inlc(void);

void  vdw_d2_stress(Lattice &, std::vector<ION> &atoms, double *stress_d2);
std::complex<double> DnmTransform(int ih, int jh, int is1, int is2, double *Ia, SPECIES &sp);
static void print_stress(char *w, double *stress_term);

template Stress<double>::~Stress(void);
template Stress<std::complex<double>>::~Stress(void);
template <class T> Stress<T>::~Stress(void)
{
}

template Stress<double>::Stress(Kpoint<double> **Kpin, Lattice &L, BaseGrid &BG, Pw &pwaves, 
        std::vector<ION> &atoms, std::vector<SPECIES> &species, double Exc, double *vxc, double *rho, double *rhocore, double *veff);

template Stress<std::complex<double>>::Stress(Kpoint<std::complex<double>> **Kpin, Lattice &L, BaseGrid &BG, Pw &pwaves, 
        std::vector<ION> &atoms, std::vector<SPECIES> &species, double Exc, double *vxc, double *rho, double *rhocore, double *veff);
template <class T> Stress<T>::Stress(Kpoint<T> **Kpin, Lattice &L, BaseGrid &BG, Pw &pwaves, 
        std::vector<ION> &atoms, std::vector<SPECIES> &species, double Exc, double *vxc, double *rho, double *rhocore, double *veff)
{

    int inlc = get_inlc();
    if(inlc !=0) 
    {
        throw RmgFatalException() << "stress does not work for vdw now" << __FILE__ << " at line " << __LINE__ << "\n";
    }

    RmgTimer *RT1 = new RmgTimer("2-Stress");
    RmgTimer *RT2;
    for(int i = 0; i < 9; i++) stress_tensor[i] = 0.0;
    RT2 = new RmgTimer("2-Stress: kinetic");
    Kinetic_term(Kpin, BG, L);
    delete RT2;
    RT2 = new RmgTimer("2-Stress: Loc");
    Local_term(atoms, species, rho, pwaves);
    delete RT2;
    RT2 = new RmgTimer("2-Stress: Non-loc");
    NonLocal_term(Kpin, atoms, species);
    delete RT2;
    if(!ct.norm_conserving_pp)  
    {
        RT2 = new RmgTimer("2-Stress: Non-loc");
        NonLocalQfunc_term(Kpin, atoms, species, veff, vxc);
        delete RT2;

    }
    RT2 = new RmgTimer("2-Stress: Hartree");
    Hartree_term(rho, pwaves);
    delete RT2;
    RT2 = new RmgTimer("2-Stress: XC Local");
    Exc_term(Exc, vxc, rho);
    delete RT2;

    RT2 = new RmgTimer("2-Stress: XC gradient");
    Exc_gradcorr(Exc, vxc, rho, rhocore);
    delete RT2;

    RT2 = new RmgTimer("2-Stress: XC nlcc");
    Exc_Nlcc(vxc, rhocore);
    delete RT2;
    RT2 = new RmgTimer("2-Stress: Ewald");
    Ewald_term(atoms, species, L, pwaves);
    delete RT2;

    if(ct.vdw_corr == DFT_D2)
    {
        RT2 = new RmgTimer("2-Stress: vdw D2 correction");
        double stress_d2[9];
        vdw_d2_stress(L, atoms, stress_d2);
        for(int i = 0; i < 9; i++) stress_tensor[i] += stress_d2[i];

        delete RT2;
    }

    if(ct.vdw_corr == DFT_D3)
        for(int i = 0; i < 9; i++) stress_tensor[i] += ct.stress_vdw[i];

    if(!ct.is_gamma && Rmg_Symm) Rmg_Symm->symmetrize_tensor(stress_tensor);
    delete RT1;
    print_stress("total ", stress_tensor);

    for(int i = 0; i < 9; i++) Rmg_L.stress_tensor[i] = stress_tensor[i];
    double zero(0.0);
    int ithree = 3;
    double  a[9]; // b is the reciprocal vector without 2PI, a^-1
    for (int i = 0; i < 3; i++)
    {
        // b[0 * 3 + i] = Rmg_L.b0[i];
        // b[1 * 3 + i] = Rmg_L.b1[i];
        // b[2 * 3 + i] = Rmg_L.b2[i];
        a[0 * 3 + i] = Rmg_L.a0[i];
        a[1 * 3 + i] = Rmg_L.a1[i];
        a[2 * 3 + i] = Rmg_L.a2[i];
    }


    double alpha = Rmg_L.omega;
    dgemm("N","N", &ithree, &ithree, &ithree, &alpha, a, &ithree, stress_tensor, &ithree, &zero, Rmg_L.cell_force, &ithree);

}

template void Stress<double>::Kinetic_term(Kpoint<double> **Kpin, BaseGrid &BG, Lattice &L);
template void Stress<std::complex<double>>::Kinetic_term(Kpoint<std::complex<double>> **Kpin, BaseGrid &BG, Lattice &L);
template <class T> void Stress<T>::Kinetic_term(Kpoint<T> **Kpin, BaseGrid &BG, Lattice &L)
{
    int PX0_GRID = Rmg_G->get_PX0_GRID(1);
    int PY0_GRID = Rmg_G->get_PY0_GRID(1);
    int PZ0_GRID = Rmg_G->get_PZ0_GRID(1);

    int pbasis = PX0_GRID * PY0_GRID * PZ0_GRID;
    int pbasis_noncol = pbasis * ct.noncoll_factor;
    T *grad_psi = (T *)GpuMallocManaged(3*pbasis_noncol*sizeof(T));
    T *psi_x = grad_psi;
    T *psi_y = psi_x + pbasis_noncol;
    T *psi_z = psi_x + 2*pbasis_noncol;

    double vel = L.get_omega() / ((double)(BG.get_NX_GRID(1) * BG.get_NY_GRID(1) * BG.get_NZ_GRID(1)));
    T alpha;
    T one(1.0);
    T stress_tensor_T[9];
    double stress_tensor_R[9];


    for(int i = 0; i < 9; i++) stress_tensor_T[i] = 0.0;

    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++)
    {
        Kpoint<T> *kptr = Kpin[kpt];
        for(int st = 0; st < kptr->nstates; st++)
        {
            if (std::abs(kptr->Kstates[st].occupation[0]) < 1.0e-10) break;
            ApplyGradient(kptr->Kstates[st].psi, psi_x, psi_y, psi_z, ct.force_grad_order, "Coarse");
            if(ct.noncoll)
                ApplyGradient(kptr->Kstates[st].psi+pbasis, psi_x+pbasis, psi_y+pbasis, psi_z+pbasis, ct.force_grad_order, "Coarse");

            if(!ct.is_gamma)
            {
                std::complex<double> I_t(0.0, 1.0);
                std::complex<double> *psi_C, *psi_xC, *psi_yC, *psi_zC;
                psi_C = (std::complex<double> *) kptr->Kstates[st].psi;
                psi_xC = (std::complex<double> *) psi_x;
                psi_yC = (std::complex<double> *) psi_y;
                psi_zC = (std::complex<double> *) psi_z;
                for(int i = 0; i < pbasis_noncol; i++)
                {
                    psi_xC[i] += I_t *  kptr->kp.kvec[0] * psi_C[i];
                    psi_yC[i] += I_t *  kptr->kp.kvec[1] * psi_C[i];
                    psi_zC[i] += I_t *  kptr->kp.kvec[2] * psi_C[i];
                }
            }


            alpha = vel * kptr->Kstates[st].occupation[0] * kptr->kp.kweight;
            RmgGemm("C", "N", 3, 3, pbasis_noncol, alpha, grad_psi, pbasis_noncol,
                    grad_psi, pbasis_noncol, one, stress_tensor_T, 3);

        }
    }

    for(int i = 0; i < 9; i++) stress_tensor_R[i] = std::real(stress_tensor_T[i])/L.omega;
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_R, 9, MPI_DOUBLE, MPI_SUM, pct.img_comm);
    if(!ct.is_gamma && Rmg_Symm) Rmg_Symm->symmetrize_tensor(stress_tensor_R);
    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_R[i];

    if(ct.verbose) print_stress("Kinetic term", stress_tensor_R);
    GpuFreeManaged(grad_psi);
}

template void Stress<double>::Hartree_term(double *rho, Pw &pwaves);
template void Stress<std::complex<double>>::Hartree_term(double *rho, Pw &pwaves);
template <class T> void Stress<T>::Hartree_term(double *rho, Pw &pwaves)
{
    int pbasis = pwaves.pbasis;
    std::complex<double> ZERO_t(0.0, 0.0);
    std::complex<double> *crho = new std::complex<double>[pbasis];


    for(int i = 0;i < pbasis;i++) crho[i] = std::complex<double>(rho[i], 0.0);
    pwaves.FftForward(crho, crho);
    for(int i = 0;i < pbasis;i++) crho[i] /=(double)pwaves.global_basis;

    double stress_tensor_h[9];
    for(int i = 0; i < 9; i ++) stress_tensor_h[i] = 0.0;
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;
    for(int ig=0;ig < pbasis;ig++) {
        if((pwaves.gmags[ig] > 1.0e-6) && pwaves.gmask[ig])
        {
            double gsqr = pwaves.gmags[ig] *tpiba2;
            double rhog2 = std::norm(crho[ig])/gsqr;
            double g[3];
            g[0] = pwaves.g[ig].a[0] * tpiba;
            g[1] = pwaves.g[ig].a[1] * tpiba;
            g[2] = pwaves.g[ig].a[2] * tpiba;
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    stress_tensor_h[i*3+j] += 2.0 * rhog2* g[i] * g[j] /gsqr;
                    if(i == j) 
                        stress_tensor_h[i*3+j] -= rhog2;
                }
            }


        }
    }

    for(int i = 0; i < 9; i++) stress_tensor_h[i] *= -0.5 * (4.0 * PI);
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_h, 9, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    if(ct.verbose) print_stress("Hartree term", stress_tensor_h);
    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_h[i];
    delete [] crho;
}

template void Stress<double>::Exc_term(double Exc, double *vxc, double *rho);
template void Stress<std::complex<double>>::Exc_term(double Exc, double *vxc, double *rho);
template <class T> void Stress<T>::Exc_term(double Exc, double *vxc, double *rho)
{
    double stress_tensor_x[9];
    for(int i = 0; i < 9; i++) stress_tensor_x[i] = 0.0;
    for(int i = 0; i < 3; i++) stress_tensor_x[i * 3 + i ] = -(ct.XC - ct.vtxc)/Rmg_L.omega;

    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_x[i];

    if(ct.verbose) print_stress("XC term", stress_tensor_x);
}

template void Stress<double>::Exc_gradcorr(double Exc, double *vxc, double *rho, double *rhocore);
template void Stress<std::complex<double>>::Exc_gradcorr(double Exc, double *vxc, double *rho, double *rhocore);
template <class T> void Stress<T>::Exc_gradcorr(double Exc, double *vxc, double *rho, double *rhocore)
{
    double stress_tensor_xcgrad[9];
    for(int i = 0; i < 9; i++) stress_tensor_xcgrad[i] = 0.0;
    Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
    if(!F->dft_is_gradient_rmg() ) return;

    int grid_ratio = Rmg_G->default_FG_RATIO;
    int pbasis = Rmg_G->get_P0_BASIS(grid_ratio);

    int fac = 1;
    if(ct.nspin == 2 ) fac = 2;
    double *rho_grad = new double[fac*3*pbasis];
    double *rho_gx = rho_grad;
    double *rho_gy = rho_grad + pbasis ;
    double *rho_gz = rho_grad + 2 * pbasis;

    double vel = Rmg_L.get_omega() / ((double)(Rmg_G->get_NX_GRID(grid_ratio) * 
                Rmg_G->get_NY_GRID(grid_ratio) * Rmg_G->get_NZ_GRID(grid_ratio)));

    double alpha = 1.0;
    if(ct.nspin == 2) alpha = 0.5;
    double malpha = -alpha;
    int ione = 1;

    F->v_xc(rho, rhocore, ct.XC, ct.vtxc, vxc, ct.nspin );

    // get gradient of rho+rhocore
    daxpy (&pbasis, &alpha, rhocore, &ione, rho, &ione);
    ApplyGradient(rho, rho_gx, rho_gy, rho_gz, ct.force_grad_order, "Fine");
    daxpy (&pbasis, &malpha, rhocore, &ione, rho, &ione);

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            for(int idx = 0; idx < pbasis; idx++)
                stress_tensor_xcgrad[i*3+j] += rho_grad[i*pbasis + idx] * rho_grad[j*pbasis + idx] * F->vxc2[idx]; 
        }

    if(ct.nspin == 2)
    {
        //opposite spin's rho is stored in &rho[pbasis];
        daxpy (&pbasis, &alpha, rhocore, &ione, &rho[pbasis], &ione);
        double *rho_gx1 = &rho_grad[3*pbasis];
        double *rho_gy1 = &rho_grad[4*pbasis];
        double *rho_gz1 = &rho_grad[5*pbasis];
        ApplyGradient(&rho[pbasis], rho_gx1, rho_gy1, rho_gz1, ct.force_grad_order, "Fine");
        daxpy (&pbasis, &malpha, rhocore, &ione, &rho[pbasis], &ione);


        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
            {
                for(int idx = 0; idx < pbasis; idx++)
                    stress_tensor_xcgrad[i*3+j] += rho_grad[i*pbasis + idx] * rho_grad[(3 + j) *pbasis + idx] * F->vxc2[idx]; 
            }

    }
    delete F;

    for(int i = 0; i < 9; i++) stress_tensor_xcgrad[i] *= vel /Rmg_L.omega;
    //  sum over grid communicator and spin communicator  but no kpoint communicator
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_xcgrad, 9, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_xcgrad, 9, MPI_DOUBLE, MPI_SUM, pct.spin_comm);
    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_xcgrad[i];

    if(ct.verbose) print_stress("XC gradcorr", stress_tensor_xcgrad);
}

template void Stress<double>::Local_term(std::vector<ION> &atoms, 
        std::vector<SPECIES> &speices, double *rho, Pw &pwaves);
template void Stress<std::complex<double>>::Local_term(std::vector<ION> &atoms, 
        std::vector<SPECIES> &speices, double *rho, Pw &pwaves);
template <class T> void Stress<T>::Local_term(std::vector<ION> &atoms, 
        std::vector<SPECIES> &speices, double *rho, Pw &pwaves)
{

    int pbasis = pwaves.pbasis;
    std::complex<double> ZERO_t(0.0, 0.0);
    std::complex<double> *crho = new std::complex<double>[pbasis];

    for(int i = 0;i < pbasis;i++) crho[i] = std::complex<double>(rho[i], 0.0);
    pwaves.FftForward(crho, crho);
    for(int i = 0;i < pbasis;i++) crho[i] /=(double)pwaves.global_basis;


    double stress_tensor_loc[9];
    for(int i = 0; i < 9; i++) stress_tensor_loc[i] = 0.0;

    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;
    double dvloc_g2, vloc_g;
    for (int isp = 0; isp < ct.num_species; isp++)
    {
        SPECIES *sp = &Species[isp];

        for(int ig=0;ig < pbasis;ig++) 
        {
            if(!pwaves.gmask[ig]) continue;
            //if((pwaves.gmags[ig] < 1.0e-6) || !pwaves.gmask[ig]) continue;
            double gsqr = pwaves.gmags[ig] *tpiba2;
            double gval = std::sqrt(gsqr);
            double g[3];
            g[0] = pwaves.g[ig].a[0] * tpiba;
            g[1] = pwaves.g[ig].a[1] * tpiba;
            g[2] = pwaves.g[ig].a[2] * tpiba;

            vloc_g = AtomicInterpolateInline_Ggrid(sp->localpp_g, gval);
            dvloc_g2 = AtomicInterpolateInline_Ggrid(sp->der_localpp_g, gval);

            std::complex<double> S;
            double gr;
            S = 0.0;
            for (int i = 0; i < ct.num_ions; i++)
            {

                ION *iptr1 = &Atoms[i];
                if (iptr1->species != isp) continue;

                gr = iptr1->crds[0] * g[0] + iptr1->crds[1] * g[1] + iptr1->crds[2] * g[2];
                S +=  std::exp(std::complex<double>(0.0, -gr));
            }

            double sg = std::real(S * std::conj(crho[ig]) );
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    stress_tensor_loc[i*3+j] += 2.0 * g[i] * g[j] * sg * dvloc_g2;
                    if(i == j) 
                        stress_tensor_loc[i*3+j] += sg * vloc_g;
                }
            }

        }
    }

    for(int i = 0; i < 9; i++) stress_tensor_loc[i] /= Rmg_L.omega;
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_loc, 9, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_loc[i];
    if(ct.verbose) print_stress("Loc term", stress_tensor_loc);
}

template void Stress<double>::NonLocal_term(Kpoint<double> **Kpin,
        std::vector<ION> &atoms, std::vector<SPECIES> &speices);
template void Stress<std::complex<double>>::NonLocal_term(Kpoint<std::complex<double>> **Kpin,
        std::vector<ION> &atoms, std::vector<SPECIES> &speices);
template <class T> void Stress<T>::NonLocal_term(Kpoint<T> **Kptr, 
        std::vector<ION> &atoms, std::vector<SPECIES> &speices)
{


    T *psi, *psi_x, *psi_y, *psi_z;
    double stress_tensor_nl[9];
    for(int i = 0; i < 9; i++) stress_tensor_nl[i] = 0.0;
    int num_occupied;
    std::complex<double> I_t(0.0, 1.0);

    int num_nonloc_ions = Kptr[0]->BetaProjector->get_num_nonloc_ions();
    int num_owned_ions = Kptr[0]->BetaProjector->get_num_owned_ions();
    int *nonloc_ions_list = Kptr[0]->BetaProjector->get_nonloc_ions_list();
    int *owned_ions_list = Kptr[0]->BetaProjector->get_owned_ions_list();


    RmgTimer *RT1;

    size_t size = num_nonloc_ions * ct.state_block_size * ct.max_nl * ct.noncoll_factor; 
    size += 1;
    //  sint_der:  leading dimension is num_nonloc_ions * ct.max_nl, 
    T *sint_der = (T *)GpuMallocManaged(size * sizeof(T));

    int num_proj = num_nonloc_ions * ct.max_nl;
    int num_proj_noncoll = num_proj * ct.noncoll_factor;
    T *proj_mat = (T *)GpuMallocManaged(num_proj_noncoll * num_proj_noncoll * sizeof(T));

    //  proj_mat_q: for US pseudopotenital only = sum_i  <beta_n *r[] |partial_ psi_i> eig[i] <psi_i|beta_n>
    T *proj_mat_q = (T *)GpuMallocManaged(num_proj_noncoll * num_proj_noncoll * sizeof(T));

    //  determine the number of occupied states for all kpoints.

    num_occupied = 0;
    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        for(int st = 0; st < ct.num_states; st++)
        {
            if(abs(Kptr[kpt]->Kstates[st].occupation[0]) < 1.0e-10) 
            {
                num_occupied = std::max(num_occupied, st);
                break;
            }
        }
    }


    int num_state_block = (num_occupied +ct.state_block_size -1) / ct.state_block_size;
    int *state_start = new int[num_state_block];
    int *state_end = new int[num_state_block];
    for(int ib = 0; ib < num_state_block; ib++)
    {
        state_start[ib] = ib * ct.state_block_size;
        state_end[ib] = (ib+1) * ct.state_block_size;
        if(state_end[ib] > num_occupied) state_end[ib] = num_occupied;
    }

    int pbasis = Kptr[0]->pbasis;
    int pbasis_noncoll = pbasis * ct.noncoll_factor;


    if(ct.alloc_states < ct.num_states + 3 * ct.state_block_size)     
    {
        printf("\n allco_states %d should be larger than ct.num_states %d + 3* ct.state_block_size, %d", ct.alloc_states, ct.num_states,
                ct.state_block_size);
        exit(0);
    }

    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {

        Kpoint<T> *kptr = Kptr[kpt];
        for(int ib = 0; ib < num_state_block; ib++)
        {
            for(int st = state_start[ib]; st < state_end[ib]; st++)
            {
                psi = Kptr[kpt]->Kstates[st].psi;
                psi_x = Kptr[kpt]->Kstates[ct.num_states].psi + (st-state_start[ib]) * pbasis_noncoll;
                psi_y = psi_x + ct.state_block_size*pbasis_noncoll;
                psi_z = psi_x +2* ct.state_block_size*pbasis_noncoll;
                ApplyGradient(psi, psi_x, psi_y, psi_z, ct.force_grad_order, "Coarse");
                if(ct.noncoll)
                    ApplyGradient(psi+pbasis, psi_x+pbasis, psi_y+pbasis, psi_z+pbasis, ct.force_grad_order, "Coarse");

                if(!ct.is_gamma)
                {
                    std::complex<double> *psi_C, *psi_xC, *psi_yC, *psi_zC;
                    psi_C = (std::complex<double> *) psi;
                    psi_xC = (std::complex<double> *) psi_x;
                    psi_yC = (std::complex<double> *) psi_y;
                    psi_zC = (std::complex<double> *) psi_z;
                    for(int i = 0; i < pbasis_noncoll; i++) 
                    {
                        psi_xC[i] += I_t *  Kptr[kpt]->kp.kvec[0] * psi_C[i];
                        psi_yC[i] += I_t *  Kptr[kpt]->kp.kvec[1] * psi_C[i];
                        psi_zC[i] += I_t *  Kptr[kpt]->kp.kvec[2] * psi_C[i];
                    }
                }

            }


            int num_state_thisblock = state_end[ib] - state_start[ib];

            for(int id1 = 0; id1 < 3; id1++)
                for(int id2 = 0; id2 < 3; id2++)
                {


                    RT1 = new RmgTimer("2-Stress: Non-loc: betaxpsi");

                    // kptr's wavefunction storage store the gradient of psi starting at first_state
                    int first_state = ct.num_states + id1 * ct.state_block_size;

                    // nlweight points to beta * x, beta * y, and beta * z for id2 = 0, 1, 2
                    T *nlweight = &kptr->nl_weight[ (id2 + 1) *kptr->nl_weight_size];

                    kptr->BetaProjector->project(kptr, sint_der, first_state *ct.noncoll_factor, num_state_thisblock *ct.noncoll_factor, nlweight);
                    delete RT1;

                    RT1 = new RmgTimer("2-Stress: Non-loc: <beta-psi>f(st) < psi-beta> ");

                    T *sint = &kptr->newsint_local[state_start[ib] * num_proj_noncoll];

                    for(int st = state_start[ib]; st < state_end[ib]; st++)
                    {
                        double t1 = kptr->Kstates[st].occupation[0] * kptr->kp.kweight;
                        for (int iproj = 0; iproj < num_proj_noncoll; iproj++)
                        {
                            sint_der[ (st-state_start[ib]) * num_proj_noncoll + iproj] *= 2.0*t1;
                            if(id1 == id2 && 1)
                            {
                                sint_der[ (st-state_start[ib]) * num_proj_noncoll + iproj] += 
                                    sint[ (st-state_start[ib]) * num_proj_noncoll + iproj] * t1; 
                            }
                        }
                    }

                    T one(1.0);
                    T zero(0.0);

                    RmgGemm("N", "C", num_proj_noncoll, num_proj_noncoll, num_state_thisblock, one, sint_der, num_proj_noncoll,
                            sint, num_proj_noncoll, zero, proj_mat, num_proj_noncoll); 

                    for(int st = state_start[ib]; st < state_end[ib]; st++)
                    {
                        double t1 = kptr->Kstates[st].eig[0];
                        for (int iproj = 0; iproj < num_proj_noncoll; iproj++)
                        {
                            sint_der[ (st-state_start[ib]) * num_proj_noncoll + iproj] *= t1;
                        }
                    }
                    RmgGemm("N", "C", num_proj_noncoll, num_proj_noncoll, num_state_thisblock, one, sint_der, num_proj_noncoll,
                            sint, num_proj_noncoll, zero, proj_mat_q, num_proj_noncoll); 

                    delete RT1;

                    RT1 = new RmgTimer("2-Stress: Non-loc: dnmI mat ");
                    int nion = -1;
                    for (int ion = 0; ion < num_owned_ions; ion++)
                    {
                        /*Global index of owned ion*/
                        int gion = owned_ions_list[ion];

                        /* Figure out index of owned ion in nonloc_ions_list array, store it in nion*/
                        do {

                            nion++;
                            if (nion >= num_nonloc_ions)
                            {
                                printf("\n Could not find matching entry in nonloc_ions_list for owned ion %d", gion);
                                rmg_error_handler(__FILE__, __LINE__, "Could not find matching entry in nonloc_ions_list for owned ion ");
                            }

                        } while (nonloc_ions_list[nion] != gion);

                        ION *iptr = &Atoms[gion];

                        int nh = Species[iptr->species].nh;
                        if(ct.noncoll)
                        {
                            std::complex<double> *dnmI_so = iptr->dnmI_so;
                            std::complex<double> *qnmI_so = iptr->qqq_so;

                            for(int is1 = 0; is1 < 2; is1++)
                                for(int is2 = 0; is2 < 2; is2++)
                                    for(int n = 0; n <nh; n++)
                                        for(int m = 0; m <nh; m++)
                                        {
                                            int idx1 = (is1 * 2 + is2) * nh * nh + n * nh + m;
                                            int ng = nion * ct.max_nl  + n  + is1 *num_proj;
                                            int mg = nion * ct.max_nl  + m  + is2 *num_proj;

                                            stress_tensor_nl[id1 * 3 + id2] += 
                                                std::real(dnmI_so[idx1] * proj_mat[ng * num_proj_noncoll + mg])
                                               -std::real(qnmI_so[idx1] * proj_mat_q[ng * num_proj_noncoll + mg]);
                                        }
                        }
                        else
                        {
                            double *dnmI = iptr->dnmI;
                            double *qnmI = iptr->qqq;


                            for(int n = 0; n <nh; n++)
                                for(int m = 0; m <nh; m++)
                                {
                                    int idx1 = n * nh + m;
                                    int ng = nion * ct.max_nl + n;
                                    int mg = nion * ct.max_nl + m;

                                    stress_tensor_nl[id1 * 3 + id2] += 
                                        dnmI[idx1] * std::real(proj_mat[ng * num_proj + mg]) -
                                        qnmI[idx1] * std::real(proj_mat_q[ng * num_proj + mg]);
                                }
                        }
                    } 
                    delete RT1;

                }

        }

    }

    delete [] state_end;
    delete [] state_start;
    GpuFreeManaged(proj_mat);
    GpuFreeManaged(proj_mat_q);
    GpuFreeManaged(sint_der);

    for(int i = 0; i < 9; i++) stress_tensor_nl[i] = stress_tensor_nl[i]/Rmg_L.omega;

    // img_comm includes kpoint, spin, and grid (num_owned_ions) sum
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_nl, 9, MPI_DOUBLE, MPI_SUM, pct.img_comm);

    if(!ct.is_gamma && Rmg_Symm) Rmg_Symm->symmetrize_tensor(stress_tensor_nl);
    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_nl[i];
    if(ct.verbose) print_stress("Nonlocal term", stress_tensor_nl);

}

template void Stress<double>::Ewald_term(std::vector<ION> &atoms, std::vector<SPECIES> &speices, 
        Lattice &L, Pw &pwaves);
template void Stress<std::complex<double>>::Ewald_term(std::vector<ION> &atoms, 
        std::vector<SPECIES> &speices, Lattice &L, Pw &pwaves);
template <class T> void Stress<T>::Ewald_term(std::vector<ION> &atoms, 
        std::vector<SPECIES> &speices, Lattice &L, Pw &pwaves)
{

    double r, x[3];
    ION *iptr1, *iptr2;
    int i, j;

    double t1;

    double sigma = 3.0;

    //  the ion-ion term in real space
    //  check how many unit cell to make erfc(nL/sqrt(2.0)*sigma) < tol ==1.0e-8
    //  erfc(4.06) = 9.37e-9 so we set the largest r/(sqrt(2)*sigma) = 4.06


    double rcutoff = 4.06;
    int num_cell_x = int(rcutoff * sqrt(2.0) * sigma/Rmg_L.get_xside()) + 1;
    int num_cell_y = int(rcutoff * sqrt(2.0) * sigma/Rmg_L.get_yside()) + 1;
    int num_cell_z = int(rcutoff * sqrt(2.0) * sigma/Rmg_L.get_zside()) + 1;


    // real space contribution
    double stress_tensor_rs[9];
    for(i=0; i < 9; i++) stress_tensor_rs[i] = 0.0;
    for (i = pct.gridpe; i < ct.num_ions; i+=pct.grid_npes)
    {

        iptr1 = &Atoms[i];
        double Zi = Species[iptr1->species].zvalence;

        for (j = 0; j < ct.num_ions; j++)
        {

            iptr2 = &Atoms[j];
            double Zj = Species[iptr2->species].zvalence;
            t1 = sqrt (sigma);

            for(int ix = -num_cell_x; ix<= num_cell_x; ix++)
                for(int iy = -num_cell_y; iy<= num_cell_y; iy++)
                    for(int iz = -num_cell_z; iz<= num_cell_z; iz++)
                    {
                        x[0] = iptr1->crds[0] - iptr2->crds[0] + ix * Rmg_L.a0[0] + iy * Rmg_L.a1[0] + iz * Rmg_L.a2[0];
                        x[1] = iptr1->crds[1] - iptr2->crds[1] + ix * Rmg_L.a0[1] + iy * Rmg_L.a1[1] + iz * Rmg_L.a2[1];
                        x[2] = iptr1->crds[2] - iptr2->crds[2] + ix * Rmg_L.a0[2] + iy * Rmg_L.a1[2] + iz * Rmg_L.a2[2];
                        r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

                        if(r < 1.0e-5) continue;

                        double hprime = 2.0/sqrt(PI) * std::exp(-sigma * r * r) - boost::math::erfc(t1*r)/(t1*r); 
                        double tem = 0.5 * t1 * Zi * Zj * hprime /(r*r);

                        for(int id1 = 0; id1 < 3; id1++)
                            for(int id2 = 0; id2 < 3; id2++)
                                stress_tensor_rs[id1 * 3 + id2] += tem * x[id1] * x[id2];

                    }
        }

    }

    //   reciprocal space term

    //if(pct.gridpe == 0) ii_kspace = -ct.nel*ct.nel / sigma / 4.0;
    double stress_tensor_gs[9];
    for(i=0; i < 9; i++) stress_tensor_gs[i] = 0.0;
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;
    double gsquare, k[3], gs4sigma;
    double kr;
    std::complex<double> S;

    for(int ig=0;ig < pwaves.pbasis;ig++)
    {
        if(pwaves.gmags[ig] > 1.0e-6)
        {
            gsquare = pwaves.gmags[ig] * tpiba2;
            gs4sigma = gsquare/(4.0 * sigma);
            k[0] = pwaves.g[ig].a[0] * tpiba;
            k[1] = pwaves.g[ig].a[1] * tpiba;
            k[2] = pwaves.g[ig].a[2] * tpiba;

            S = 0.0;
            for (int i = 0; i < ct.num_ions; i++)
            {

                iptr1 = &Atoms[i];
                double Zi = Species[iptr1->species].zvalence;
                kr = iptr1->crds[0] * k[0] + iptr1->crds[1] * k[1] + iptr1->crds[2] * k[2];
                S +=  Zi * std::exp(std::complex<double>(0.0, kr));
            }

            double tem = std::norm(S) * exp(-gs4sigma)/ gs4sigma; 
            double tem1 = 2.0 * tem/gsquare *(gs4sigma +1.0);
            for(int id1 = 0; id1 < 3; id1++)
            {
                for(int id2 = 0; id2 < 3; id2++)
                {
                    stress_tensor_gs[id1 * 3 + id2] += tem1 * k[id1] * k[id2];
                }
            }
            for(int id1 = 0; id1 < 3; id1++)
            {
                stress_tensor_gs[id1 * 3 + id1] -= tem;
            }
        }
    }

    for(i=0; i < 9; i++) stress_tensor_gs[i] *= PI/(2.0 * Rmg_L.omega * sigma);
    for(i=0; i < 9; i++) stress_tensor_gs[i] += stress_tensor_rs[i];

    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_gs, 9, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    for(i = 0; i < 3; i++) stress_tensor_gs[i*3+i] += PI/(2.0*Rmg_L.omega * sigma) * ct.nel * ct.nel;

    for(i=0; i < 9; i++) stress_tensor_gs[i] /= -Rmg_L.omega;

    for(i=0; i < 9; i++) stress_tensor[i] += stress_tensor_gs[i];

    if(ct.verbose) print_stress("Ewald term", stress_tensor_gs);
}


template void Stress<double>::NonLocalQfunc_term(Kpoint<double> **Kpin,
        std::vector<ION> &atoms, std::vector<SPECIES> &speices, double *veff, double *vxc);
template void Stress<std::complex<double>>::NonLocalQfunc_term(Kpoint<std::complex<double>> **Kpin,
        std::vector<ION> &atoms, std::vector<SPECIES> &speices, double *veff, double *vxc);
template <class T> void Stress<T>::NonLocalQfunc_term(Kpoint<T> **Kptr, 
        std::vector<ION> &atoms, std::vector<SPECIES> &speices, double *veff, double *vxc)
{


    double stress_tensor_nlq[9];
    for(int i = 0; i < 9; i++) stress_tensor_nlq[i] = 0.0;
    int num_occupied;
    std::complex<double> I_t(0.0, 1.0);


    int num_nonloc_ions = Kptr[0]->BetaProjector->get_num_nonloc_ions();
    int num_owned_ions = Kptr[0]->BetaProjector->get_num_owned_ions();
    int *nonloc_ions_list = Kptr[0]->BetaProjector->get_nonloc_ions_list();
    int *owned_ions_list = Kptr[0]->BetaProjector->get_owned_ions_list();



    int num_proj = num_nonloc_ions * ct.max_nl;
    int num_proj_noncoll = num_proj * ct.noncoll_factor;
    T *proj_mat = (T *)GpuMallocManaged(num_proj_noncoll * num_proj_noncoll * sizeof(T));

    //  determine the number of occupied states for all kpoints.

    num_occupied = 0;
    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        for(int st = 0; st < ct.num_states; st++)
        {
            if(abs(Kptr[kpt]->Kstates[st].occupation[0]) < 1.0e-10) 
            {
                num_occupied = std::max(num_occupied, st);
                break;
            }
        }
    }


    size_t size = num_nonloc_ions * ct.max_nl * num_occupied * ct.noncoll_factor; 
    size += 1;
    //  sint_der:  leading dimension is num_nonloc_ions * ct.max_nl, 
    T *sint_der = (T *)GpuMallocManaged(size * sizeof(T));


    // proj_mat = Sum_stm kpt (kpweigth *  <beta_n|psi_i> occ[i] <psi_i|beta_m>) eq 27 in PRB 61, 8433
    for(int i = 0; i < num_proj_noncoll * num_proj_noncoll; i++) proj_mat[i] = 0.0;
    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {

        Kpoint<T> *kptr = Kptr[kpt];
        T *sint = kptr->newsint_local;
        for(int st = 0; st < num_occupied; st++)
        {
            double t1 = kptr->Kstates[st].occupation[0] * kptr->kp.kweight;
            for (int iproj = 0; iproj < num_proj_noncoll; iproj++)
            {
                sint_der[ st * num_proj_noncoll + iproj] = t1 * sint[ st * num_proj_noncoll + iproj];

            }
        }

        T one(1.0);

        RmgGemm("N", "C", num_proj_noncoll, num_proj_noncoll, num_occupied, one, sint_der, num_proj_noncoll,
                sint, num_proj_noncoll, one, proj_mat, num_proj_noncoll); 
    }


    int idx, i, j, ion;
    int nh, ncount, icount;
    double *sum;
    int *ivec, sum_dim;
    ION *iptr;
    SPECIES *sp;

    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    /*Count the number of elements in sum array */
    int nh_max = 0;
    for(int sp = 0; sp <ct.num_species; sp++)
        nh_max = std::max(nh_max, Species[sp].nh);
    int num_q_max =  ((nh_max+1) * nh_max)/2;
    sum_dim = ct.num_ions * num_q_max;
    int sum_dim_noncoll = sum_dim * ct.noncoll_factor * ct.noncoll_factor;

    sum = new double[sum_dim_noncoll];

    double *veff_grad = new double[3*FP0_BASIS];
    double *veff_gx = veff_grad;
    double *veff_gy = veff_grad + FP0_BASIS;
    double *veff_gz = veff_grad + 2*FP0_BASIS;
    double *veff_gxyz;

    double *vxc_grad = NULL;
    double *vxc_gx = NULL;
    double *vxc_gy = NULL;
    double *vxc_gz = NULL;
    double *vxc_gxyz = NULL;

    if(ct.noncoll)
    {
        vxc_grad = new double[3*3*FP0_BASIS];
        vxc_gx = vxc_grad;
        vxc_gy = vxc_grad + FP0_BASIS * 3;
        vxc_gz = vxc_grad + FP0_BASIS * 6;
        ApplyGradient(&vxc[FP0_BASIS], vxc_gx, vxc_gy, vxc_gz, ct.force_grad_order, "Fine");
        ApplyGradient(&vxc[2*FP0_BASIS], vxc_gx+FP0_BASIS, vxc_gy+FP0_BASIS, vxc_gz+FP0_BASIS, ct.force_grad_order, "Fine");
        ApplyGradient(&vxc[3*FP0_BASIS], vxc_gx+2*FP0_BASIS, vxc_gy+2*FP0_BASIS, vxc_gz+2*FP0_BASIS, ct.force_grad_order, "Fine");
    }

    ApplyGradient(veff, veff_gx, veff_gy, veff_gz, ct.force_grad_order, "Fine");
    for(int id1 = 0; id1 < 3; id1++)
        for(int id2 = 0; id2 < 3; id2++)
        {

            veff_gxyz = veff_grad + id1 * FP0_BASIS;
            if(ct.noncoll) vxc_gxyz = vxc_grad + id1 * 3 * FP0_BASIS;
            for(int i = 0; i < sum_dim_noncoll; i++) sum[i] = 0.0;

            for (ion = 0; ion < ct.num_ions; ion++)
            {
                iptr = &Atoms[ion];
                sp = &Species[iptr->species];

                ivec = Atoms[ion].Qindex.data();
                nh = sp->nh;
                ncount = Atoms[ion].Qindex.size();

                idx = 0;
                for (i = 0; i < nh; i++)
                {
                    for (j = i; j < nh; j++)
                    {
                        if (ncount)
                        {
                            for (icount = 0; icount < ncount; icount++)
                            {
                                sum[ion * num_q_max + idx] += Atoms[ion].augfunc_xyz[id2][icount + idx * ncount] * veff_gxyz[ivec[icount]];
                            }

                            if(ct.noncoll)
                            {
                                for (icount = 0; icount < ncount; icount++)
                                {
                                    sum[ion * num_q_max + idx + sum_dim] += 
                                        Atoms[ion].augfunc_xyz[id2][icount + idx * ncount] * vxc_gxyz[ivec[icount] + 0 * FP0_BASIS];
                                    sum[ion * num_q_max + idx + 2*sum_dim] += 
                                        Atoms[ion].augfunc_xyz[id2][icount + idx * ncount] * vxc_gxyz[ivec[icount] + 1 * FP0_BASIS];
                                    sum[ion * num_q_max + idx + 3*sum_dim] += 
                                        Atoms[ion].augfunc_xyz[id2][icount + idx * ncount] * vxc_gxyz[ivec[icount] + 2 * FP0_BASIS];
                                }
                            }
                            // switching the derivative from Q to veff introduce Q * Veff term 
                            // cancelling out the original Q * Veff term
                            if(id1 == id2 && 0 )  
                                for (icount = 0; icount < ncount; icount++)
                                {
                                    sum[ion * num_q_max + idx] += Atoms[ion].augfunc[icount + idx * ncount] * veff[ivec[icount]];
                                }

                        }               /*end if (ncount) */


                        idx++;
                    }                   /*end for (j = i; j < nh; j++) */
                }                       /*end for (i = 0; i < nh; i++) */

            }                           /*end for (ion = 0; ion < ct.num_ions; ion++) */

            for(int idx = 0; idx < sum_dim_noncoll; idx++)
                sum[idx] *= get_vel_f();
            global_sums (sum, &sum_dim_noncoll, pct.grid_comm);

            int nion = -1;
            double *Ia = new double[ct.max_nl * ct.max_nl * 4];
            std::complex<double> dnmI_so;
            for (int ion = 0; ion < num_owned_ions; ion++)
            {
                /*Global index of owned ion*/
                int gion = owned_ions_list[ion];

                /* Figure out index of owned ion in nonloc_ions_list array, store it in nion*/
                do {

                    nion++;
                    if (nion >= num_nonloc_ions)
                    {
                        printf("\n Could not find matching entry in nonloc_ions_list for owned ion %d", gion);
                        rmg_error_handler(__FILE__, __LINE__, "Could not find matching entry in nonloc_ions_list for owned ion ");
                    }

                } while (nonloc_ions_list[nion] != gion);

                ION *iptr = &Atoms[gion];
                sp = &Species[iptr->species];

                int nh = Species[iptr->species].nh;

                if(!ct.noncoll)
                {
                    int idx1 = 0;
                    for(int n = 0; n <nh; n++)
                        for(int m = n; m <nh; m++)
                        {
                            int ng = nion * ct.max_nl + n;
                            int mg = nion * ct.max_nl + m;

                            double fac = 1.0;
                            if(n != m) fac = 2.0;
                            stress_tensor_nlq[id1 * 3 + id2] += fac* sum[gion * num_q_max + idx1] * std::real(proj_mat[ng * num_proj + mg]);
                            idx1++;
                        }
                }
                else
                {


                    int idx1 = 0;
                    for (int ih = 0; ih < nh; ih++)
                        for (int jh = ih; jh < nh; jh++)
                        {
                            for(int is = 0; is < 4; is++)
                            {
                                Ia[ih * nh + jh + is *nh*nh] = sum[is * sum_dim + gion * num_q_max + idx1];
                                Ia[jh * nh + ih + is *nh*nh] = sum[is * sum_dim + gion * num_q_max + idx1];
                            }
                            idx1++;
                        }

                    for(int is1 = 0; is1 <2; is1++)
                        for(int is2 = 0; is2 <2; is2++)
                        {

                            for (int ih = 0; ih < nh; ih++)
                            {
                                for (int jh = 0; jh < nh; jh++)
                                {
                                    int ng = nion * ct.max_nl + ih + is1 * num_proj;
                                    int mg = nion * ct.max_nl + jh + is2 * num_proj;
                                    dnmI_so = DnmTransform(ih,jh,is1,is2,Ia,*sp);
                                    stress_tensor_nlq[id1 * 3 + id2] += std::real(dnmI_so * proj_mat[ng * num_proj_noncoll + mg]);
                                }
                            }           
                        }              

                }

            }
            delete [] Ia;
        }

    for(int i = 0; i < 9; i++) stress_tensor_nlq[i] = stress_tensor_nlq[i]/Rmg_L.omega;

    // img_comm includes spin, and grid (num_owned_ions) sum
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_nlq, 9, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_nlq, 9, MPI_DOUBLE, MPI_SUM, pct.spin_comm);

    if(!ct.is_gamma && Rmg_Symm) Rmg_Symm->symmetrize_tensor(stress_tensor_nlq);
    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_nlq[i];
    if(ct.verbose) print_stress("NonlocalQfunc term", stress_tensor_nlq);
    delete [] sum;
    delete [] veff_grad;
    GpuFreeManaged(sint_der);
    GpuFreeManaged(proj_mat);


}

template void Stress<double>::Exc_Nlcc(double *vxc, double *rhocore);
template void Stress<std::complex<double>>::Exc_Nlcc(double *vxc, double *rhocore);
template <class T> void Stress<T>::Exc_Nlcc(double *vxc, double *rhocore)
{
    double stress_tensor_xcc[9];
    for(int i = 0; i < 9; i++) stress_tensor_xcc[i] = 0.0;

    int grid_ratio = Rmg_G->default_FG_RATIO;
    int pbasis = Rmg_G->get_P0_BASIS(grid_ratio);

    int fac = 1.0;
    if(ct.spin_flag ) fac = 2;
    double *vxc_grad = new double[fac*3*pbasis];
    double *vxc_gx = vxc_grad;
    double *vxc_gy = vxc_grad + pbasis ;
    double *vxc_gz = vxc_grad + 2 * pbasis;
    double *rhocore_stress = new double[fac*3*pbasis];

    double vel = Rmg_L.get_omega() / ((double)(Rmg_G->get_NX_GRID(grid_ratio) * 
                Rmg_G->get_NY_GRID(grid_ratio) * Rmg_G->get_NZ_GRID(grid_ratio)));

    double *dummy = NULL;
    InitLocalObject(rhocore_stress,dummy, ATOMIC_RHOCORE_STRESS, false); 
    double alpha = vel;
    // for spin-polarized case, the rhocore should be split into half+half
    if(ct.spin_flag) alpha = 0.5 * vel;
    double zero = 0.0;

    ApplyGradient(vxc, vxc_gx, vxc_gy, vxc_gz, ct.force_grad_order, "Fine");

    int ithree = 3;
    dgemm("T","N", &ithree, &ithree, &pbasis, &alpha, vxc_grad, &pbasis, rhocore_stress, &pbasis, &zero, stress_tensor_xcc, &ithree);

    double diag_term = 0.0;
    for(int idx = 0; idx < pbasis; idx++) diag_term += alpha * vxc[idx] * rhocore[idx];

    for(int i = 0; i < 3; i++) stress_tensor_xcc[i*3+i] += diag_term;

    for(int i = 0; i < 9; i++) stress_tensor_xcc[i] *= 1.0 /Rmg_L.omega;
    //  sum over grid communicator and spin communicator  but no kpoint communicator
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_xcc, 9, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_xcc, 9, MPI_DOUBLE, MPI_SUM, pct.spin_comm);
    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_xcc[i];

    if(ct.verbose) print_stress("XC Nlcc", stress_tensor_xcc);
    delete [] vxc_grad;
    delete [] rhocore_stress;
}

static void print_stress(char *w, double *stress_term)
{
    if(pct.imgpe == 0)
    {
        printf("\n stress  %s in unit of kbar", w);
        fprintf(ct.logfile, "\n stress %s  in unit of kbar", w);
        for(int i = 0; i < 3; i++)
        {
            printf("\n");
            fprintf(ct.logfile, "\n");
            for (int j = 0; j < 3; j++) printf(" %f ", stress_term[i*3+j] * Ha_Kbar);
            for (int j = 0; j < 3; j++) fprintf(ct.logfile, " %f ", stress_term[i*3+j] * Ha_Kbar);
        }
        printf("\n");
        fprintf(ct.logfile, "\n");
    }
    fflush(NULL);
}

