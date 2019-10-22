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


static void dlocal_pp(double gval, int rg_points, double *r, double *rab, double *vloc0, double *work, 
        double &dvloc_g, double &vloc_g);
static void print_stress(char *w, double *stress_term);

template Stress<double>::~Stress(void);
template Stress<std::complex<double>>::~Stress(void);
template <class T> Stress<T>::~Stress(void)
{
}

template Stress<double>::Stress(Kpoint<double> **Kpin, Lattice &L, BaseGrid &BG, Pw &pwaves, 
        std::vector<ION> &atoms, std::vector<SPECIES> &species, double Exc, double *vxc, double *rho);
template Stress<std::complex<double>>::Stress(Kpoint<std::complex<double>> **Kpin, Lattice &L, BaseGrid &BG, Pw &pwaves, 
        std::vector<ION> &atoms, std::vector<SPECIES> &species, double Exc, double *vxc, double *rho);
template <class T> Stress<T>::Stress(Kpoint<T> **Kpin, Lattice &L, BaseGrid &BG, Pw &pwaves, 
        std::vector<ION> &atoms, std::vector<SPECIES> &species, double Exc, double *vxc, double *rho)
{

    RmgTimer *RT1 = new RmgTimer("2-Stress");
    RmgTimer *RT2;
    for(int i = 0; i < 9; i++) stress_tensor[i] = 0.0;
    RT2 = new RmgTimer("2-Stress: kinetic");
    Kinetic_term(Kpin, BG, L);
    delete RT2;
    RT2 = new RmgTimer("2-Stress: Loc");
    Local_term(atoms, species, rho, pwaves);
    //Local_term1(rho, vnuc);
    delete RT2;
    RT2 = new RmgTimer("2-Stress: Non-loc");
    NonLocal_term(Kpin, atoms, species);
    delete RT2;
    RT2 = new RmgTimer("2-Stress: Hartree");
    Hartree_term(rho, pwaves);
    delete RT2;
    RT2 = new RmgTimer("2-Stress: XC");
    Exc_term(Exc, vxc, rho);
    delete RT2;
    RT2 = new RmgTimer("2-Stress: Ewald");
    Ewald_term(atoms, species, L, pwaves);
    delete RT2;
    delete RT1;
    print_stress("total stress", stress_tensor);
}

template void Stress<double>::Kinetic_term(Kpoint<double> **Kpin, BaseGrid &BG, Lattice &L);
template void Stress<std::complex<double>>::Kinetic_term(Kpoint<std::complex<double>> **Kpin, BaseGrid &BG, Lattice &L);
template <class T> void Stress<T>::Kinetic_term(Kpoint<T> **Kpin, BaseGrid &BG, Lattice &L)
{
    int PX0_GRID = Rmg_G->get_PX0_GRID(1);
    int PY0_GRID = Rmg_G->get_PY0_GRID(1);
    int PZ0_GRID = Rmg_G->get_PZ0_GRID(1);

    int pbasis = PX0_GRID * PY0_GRID * PZ0_GRID;
    T *grad_psi = (T *)GpuMallocManaged(3*pbasis*sizeof(T));
    T *psi_x = grad_psi;
    T *psi_y = psi_x + pbasis;
    T *psi_z = psi_x + 2*pbasis;

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
            //ApplyGradient(kptr->Kstates[st].psi, psi_x, psi_y, psi_z, 0, "Coarse");

            alpha = vel * kptr->Kstates[st].occupation[0];
            RmgGemm("C", "N", 3, 3, pbasis, alpha, grad_psi, pbasis,
                    grad_psi, pbasis, one, stress_tensor_T, 3);
        }
    }

    for(int i = 0; i < 9; i++) stress_tensor_R[i] = std::real(stress_tensor_T[i])/L.omega;
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_R, 9, MPI_DOUBLE, MPI_SUM, pct.img_comm);
    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_R[i];

    print_stress("Kinetic term", stress_tensor_R);
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
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_h, 9, MPI_DOUBLE, MPI_SUM, pct.img_comm);
    print_stress("Hartree term", stress_tensor_h);
    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_h[i];
    delete [] crho;
}

template void Stress<double>::Exc_term(double Exc, double *vxc, double *rho);
template void Stress<std::complex<double>>::Exc_term(double Exc, double *vxc, double *rho);
template <class T> void Stress<T>::Exc_term(double Exc, double *vxc, double *rho)
{
    int PX0_GRID = Rmg_G->get_PX0_GRID(1);
    int PY0_GRID = Rmg_G->get_PY0_GRID(1);
    int PZ0_GRID = Rmg_G->get_PZ0_GRID(1);

    int pbasis = PX0_GRID * PY0_GRID * PZ0_GRID;
    double stress_tensor_x[9];
    for(int i = 0; i < 9; i++) stress_tensor_x[i] = 0.0;
    double tem = 0.0;
    for(int i = 0; i < 3; i++) stress_tensor_x[i * 3 + i ] = -(ct.XC - ct.vtxc)/Rmg_L.omega;

    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_x[i];

    print_stress("XC term", stress_tensor_x);
}

template void Stress<double>::Local_term1(double *rho, double *vnuc);
template void Stress<std::complex<double>>::Local_term1(double *rho, double *vnuc);
template <class T> void Stress<T>::Local_term1(double *rho, double *vnuc)
{

    int grid_ratio = Rmg_G->default_FG_RATIO;
    int PX0_GRID = Rmg_G->get_PX0_GRID(grid_ratio);
    int PY0_GRID = Rmg_G->get_PY0_GRID(grid_ratio);
    int PZ0_GRID = Rmg_G->get_PZ0_GRID(grid_ratio);

    int pbasis = PX0_GRID * PY0_GRID * PZ0_GRID;
    double *grad_vnuc = new double[3 * pbasis];
    double *vnuc_x = grad_vnuc;
    double *vnuc_y = grad_vnuc + pbasis;
    double *vnuc_z = grad_vnuc + 2*pbasis;

    double vel = Rmg_L.get_omega() / ((double)(Rmg_G->get_NX_GRID(grid_ratio) * Rmg_G->get_NY_GRID(grid_ratio) * Rmg_G->get_NZ_GRID(grid_ratio)));
    double stress_tensor_loc[9];


    for(int i = 0; i < 9; i++) stress_tensor_loc[i] = 0.0;
    ApplyGradient(vnuc, vnuc_x, vnuc_y, vnuc_z, ct.force_grad_order, "Fine");

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
        double *work = new double[sp->rg_points];
        double *vloc0 = new double[sp->rg_points];

        double Zv = sp->zvalence;
        double fac = 4.0 * PI * Zv;
        double rc = sp->rc;
        // add the long range compensating charge contribution in real space
        for (int idx = 0; idx < sp->rg_points; idx++)
            vloc0[idx] = sp->vloc0[idx] + Zv * boost::math::erf (sp->r[idx] / rc) / sp->r[idx];

        for(int ig=0;ig < pbasis;ig++) 
        {
            if((pwaves.gmags[ig] < 1.0e-6) || !pwaves.gmask[ig]) continue;
            //if((pwaves.gmags[ig] < 1.0e-6) ) continue;
            double gsqr = pwaves.gmags[ig] *tpiba2;
            double gval = std::sqrt(gsqr);
            double g[3];
            g[0] = pwaves.g[ig].a[0] * tpiba;
            g[1] = pwaves.g[ig].a[1] * tpiba;
            g[2] = pwaves.g[ig].a[2] * tpiba;

            //dlocal_pp(gval, sp->rg_points, sp->r, sp->rab, vloc0, work, dvloc_g2, vloc_g);

            // substract the long range part in G space.
            //vloc_g -= fac * std::exp(-0.25 * gsqr *rc *rc)/gsqr;
            //dvloc_g2 += fac * std::exp(-0.25 * gsqr *rc *rc)/(gsqr*gsqr) *( 0.25 * gsqr * rc * rc + 1.0);
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
                S +=  std::exp(std::complex<double>(0.0, gr));
            }

            double sg = std::real(S * std::conj(crho[ig]) );
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    stress_tensor_loc[i*3+j] -= 2.0 * g[i] * g[j] * sg * dvloc_g2;
                    if(i == j) 
                        stress_tensor_loc[i*3+j] -= sg * vloc_g;
                }
            }

        }
    }

    for(int i = 0; i < 9; i++) stress_tensor_loc[i] /= Rmg_L.omega;
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_loc, 9, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_loc[i];
    print_stress("Loc term", stress_tensor_loc);
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

    int P0_BASIS = Rmg_G->get_P0_BASIS(1);

    int num_nonloc_ions = Kptr[0]->BetaProjector->get_num_nonloc_ions();
    int num_owned_ions = Kptr[0]->BetaProjector->get_num_owned_ions();
    int *nonloc_ions_list = Kptr[0]->BetaProjector->get_nonloc_ions_list();
    int *owned_ions_list = Kptr[0]->BetaProjector->get_owned_ions_list();


    RmgTimer *RT1;

    size_t size = num_nonloc_ions * ct.state_block_size * ct.max_nl; 
    size += 1;
//  sint_der:  leading dimension is num_nonloc_ions * ct.max_nl, 
    T *sint_der = (T *)GpuMallocManaged(size * sizeof(T));

    int num_proj = num_nonloc_ions * ct.max_nl;
    T *proj_mat = (T *)GpuMallocManaged(num_proj * num_proj * sizeof(T));

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
                psi_x = Kptr[kpt]->Kstates[ct.num_states].psi + (st-state_start[ib]) * pbasis;
                psi_y = psi_x + ct.state_block_size*pbasis;
                psi_z = psi_x +2* ct.state_block_size*pbasis;
                ApplyGradient(psi, psi_x, psi_y, psi_z, ct.force_grad_order, "Coarse");
                //ApplyGradient(psi, psi_x, psi_y, psi_z, 0, "Coarse");

                if(!ct.is_gamma)
                {
                    std::complex<double> *psi_C, *psi_xC, *psi_yC, *psi_zC;
                    psi_C = (std::complex<double> *) psi;
                    psi_xC = (std::complex<double> *) psi_x;
                    psi_yC = (std::complex<double> *) psi_y;
                    psi_zC = (std::complex<double> *) psi_z;
                    for(int i = 0; i < P0_BASIS; i++) 
                    {
                        psi_xC[i] += I_t *  Kptr[kpt]->kp.kvec[0] * psi_C[i];
                        psi_yC[i] += I_t *  Kptr[kpt]->kp.kvec[1] * psi_C[i];
                        psi_zC[i] += I_t *  Kptr[kpt]->kp.kvec[2] * psi_C[i];
                    }
                }


            }


            int num_state_thisblock = state_end[ib] - state_start[ib];


            for(int id1 = 0, id2 = 0; id1 < 3 && id2 < 3; id1++, id2++)
            {

                RT1 = new RmgTimer("2-Stress: Non-loc: betaxpsi");

                // kptr's wavefunction storage store the gradient of psi starting at first_state
                int first_state = ct.num_states + id1 * ct.state_block_size;

                // nlweight points to beta * x, beta * y, and beta * z for id2 = 0, 1, 2
                T *nlweight = &kptr->nl_weight[ (id2 + 1) *kptr->nl_weight_size];

                kptr->BetaProjector->project(kptr, sint_der, first_state, num_state_thisblock, nlweight);
                delete RT1;

                RT1 = new RmgTimer("2-Stress: Non-loc: <beta-psi>f(st) < psi-beta> ");

                T *sint = &kptr->newsint_local[state_start[ib] * num_proj];

                for(int st = state_start[ib]; st < state_end[ib]; st++)
                {
                    double t1 = kptr->Kstates[st].occupation[0] * kptr->kp.kweight;
                    for (int iproj = 0; iproj < num_proj; iproj++)
                    {
                        sint_der[ (st-state_start[ib]) * num_proj + iproj] *= t1;
                        if(id1 == id2  && 0)
                        {
                            sint_der[ (st-state_start[ib]) * num_proj + iproj] += 
                                sint[ (st-state_start[ib]) * num_proj + iproj] * t1; 
                        }
                    }
                }

                T one(1.0);
                T zero(0.0);

                RmgGemm("N", "C", num_proj, num_proj, num_state_thisblock, one, sint_der, num_proj,
                        sint, num_proj, zero, proj_mat, num_proj); 


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
                    double *dnmI = pct.dnmI[gion];

                    for(int n = 0, m = 0; n <nh && m <nh; n++, m++)
                    {
                        int idx1 = n * nh + m;
                        int ng = nion * ct.max_nl + n;
                        int mg = nion * ct.max_nl + m;

                        stress_tensor_nl[id1 * 3 + id2] += 
                            dnmI[idx1] * std::real(proj_mat[ng * num_proj + mg]);
                    }
                } 
                delete RT1;

            }

        }

    }

    delete [] state_end;
    delete [] state_start;

    for(int i = 0; i < 9; i++) stress_tensor_nl[i] = stress_tensor_nl[i]/Rmg_L.omega;

    // img_comm includes kpoint, spin, and grid (num_owned_ions) sum
    MPI_Allreduce(MPI_IN_PLACE, stress_tensor_nl, 9, MPI_DOUBLE, MPI_SUM, pct.img_comm);
    for(int i = 0; i < 9; i++) stress_tensor[i] += stress_tensor_nl[i];
    print_stress("Nonlocal term", stress_tensor_nl);

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

    print_stress("Ewald term", stress_tensor_gs);
}

static void dlocal_pp(double gval, int rg_points, double *r, double *rab, 
        double *vloc0, double *work, double &dvloc_g2, double &vloc_g)
{
    //  vloc0(r) = vloc(r) + Zv * erf(r/rc) /r in radial grid 
    //  vloc_g = integration of  vloc0(r) * sin(gr)/gr 
    //  dvloc_g:  d_vloc_g/d_g = integration of vloc0(r) * (gr*cos(gr) - sin(gr))/(g^2r)
    //  dvloc_g2: d_vloc_g/d_g^2 = dvloc_g/(2.0*g)
    for(int idx = 0; idx < rg_points; idx++)
    {
        double gr = gval * r[idx];
        work[idx] = vloc0[idx] * (gr * std::cos(gr) - std::sin(gr)) /(gval * gr);
    }
    dvloc_g2 = 4.0 * PI * radint1 (work, r, rab, rg_points)/(2.0 * gval);

//    for(int idx = 0; idx < rg_points; idx++)
//    {
//        double gr = gval * r[idx];
//        work[idx] = vloc0[idx] * std::sin(gr) /gr;
//    }
//    vloc_g = 4.0 * PI * radint1 (work, r, rab, rg_points)/(2.0 * gval);
}

static void print_stress(char *w, double *stress_term)
{
    if(pct.gridpe == 0)
    {
        printf("\n stress term %s", w);
        fprintf(ct.logfile, "\n stress term %s", w);
        for(int i = 0; i < 3; i++)
        {
            printf("\n");
            fprintf(ct.logfile, "\n");
            for (int j = 0; j < 3; j++) printf(" %f ", stress_term[i*3+j] * Ha_Kbar);
            for (int j = 0; j < 3; j++) fprintf(ct.logfile, " %f ", stress_term[i*3+j] * Ha_Kbar);
        }
    }
}

