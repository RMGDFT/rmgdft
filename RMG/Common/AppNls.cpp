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


#include "transition.h"
#include "const.h"
#include "BaseThread.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "Functional.h"

#include "GlobalSums.h"
#include "blas.h"

template void AppNls<double>(Kpoint<double> *, double *, double *, double *, double *, int, int);
template void AppNls<std::complex<double> >(Kpoint<std::complex<double>> *, std::complex<double> *, 
        std::complex<double> *, std::complex<double> *, std::complex<double> *, int, int);

    template <typename KpointType>
void AppNls(Kpoint<KpointType> *kpoint, KpointType *sintR, 
        KpointType *psi, KpointType *nv, KpointType *ns,
        int first_state, int num_states)
{

    // Sanity check
    if(num_states > ct.non_local_block_size)
        throw RmgFatalException() << "AppNls called with num_states > non_local_block_size in " << __FILE__ << " at line " << __LINE__ << "\n";

    KpointType *weight = kpoint->nl_weight;
#if HIP_ENABLED || CUDA_ENABLED
    weight = kpoint->nl_weight_gpu;
#endif

    //   sintR:  <beta | psi_up, psi_down>, dimensiont is numProj * 2 * num_states in noncollinear case.
    //   nv : |beta_n> Dnm <beta_m|psi_up, psi_down>, dimension 2 * pbasis
    //   ns : |psi_up, psu_down > + |beta_n> Dnm <beta_m|psi_up, psi_down>, dimension 2 * pbasis

    int P0_BASIS = kpoint->pbasis;
    int num_nonloc_ions = kpoint->BetaProjector->get_num_nonloc_ions();
    int *nonloc_ions_list = kpoint->BetaProjector->get_nonloc_ions_list();
    int num_tot_proj = kpoint->BetaProjector->get_num_tot_proj();

    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);

    char *transa = "n";

    double *dnmI;
    double *qqq;
    size_t stop = (size_t)num_states * (size_t)P0_BASIS * (size_t) ct.noncoll_factor;

    for(size_t i = 0; i < stop; i++) nv[i] = ZERO_t;
    if(num_tot_proj == 0)
    {
        bool need_ns = true;
        if(ct.norm_conserving_pp && ct.is_gamma) need_ns = false;
        if(need_ns) for(size_t idx = 0;idx < stop;idx++) ns[idx] = psi[idx];
        if(ct.xc_is_hybrid && Functional::is_exx_active())
        {
            for(size_t i = 0; i < stop; i++) nv[i] = ct.exx_fraction * kpoint->vexx[(size_t)first_state*(size_t)P0_BASIS + i];
            //AppExx(kpoint, psi, num_states, kpoint->vexx, &nv[(size_t)first_state*(size_t)P0_BASIS]);
        }
        else
        {
            for(size_t i = 0; i < stop; i++) nv[i] = ZERO_t;
        }
        return;
    }


    size_t alloc = (size_t)num_tot_proj * (size_t)num_states * ct.noncoll_factor;
    size_t M_cols = (size_t)num_tot_proj * ct.noncoll_factor;
    size_t alloc1 = (size_t)ct.max_nl * (size_t)M_cols * ct.noncoll_factor;

    KpointType *sint_compack = (KpointType *)RmgMallocHost(sizeof(KpointType) * alloc);
    KpointType *nwork = (KpointType *)RmgMallocHost(sizeof(KpointType) * alloc);
    KpointType *nwork_ion = (KpointType *)RmgMallocHost(sizeof(KpointType) * alloc);
    KpointType *M_dnm = (KpointType *)RmgMallocHost(sizeof(KpointType) * alloc1);
    KpointType *M_qqq = (KpointType *)RmgMallocHost(sizeof(KpointType) * alloc1);
    std::complex<double> *M_dnm_C = (std::complex<double> *) M_dnm;
    std::complex<double> *M_qqq_C = (std::complex<double> *) M_qqq;

    size_t sindex = first_state * ct.noncoll_factor * num_nonloc_ions * ct.max_nl;

    for (size_t i = 0; i < alloc1; i++)
    {
        M_dnm[i] = ZERO_t;
        M_qqq[i] = ZERO_t;
    }

    for (size_t i = 0; i < alloc; i++)
    {
        sint_compack[i] = 0.0;
    }

    // sintR (ct.max_nl, num_nonloc_ions, noncoll, num_states)
    // sint_compack(ct.max_nl, noncoll, num_states, num_nonloc_ions)
    for (int ion = 0; ion < num_nonloc_ions; ion++)
    {
        size_t idx = ion * ct.max_nl * ct.noncoll_factor * num_states;
        for(int st = 0; st < ct.noncoll_factor * num_states; st++)
        {
            for(int ih = 0; ih < ct.max_nl; ih++)
            {
                sint_compack[idx + st * ct.max_nl + ih] = sintR[sindex + st * num_tot_proj + ion * ct.max_nl + ih];
            }
        }
    }

    // set up M_qqq and M_dnm, this can be done outside in the
    //  M_dnm (ct.max_nl, noncoll, ct.max_nl, noncoll, num_nonloc_ions)
    //  M_qqq (ct.max_nl, noncoll, ct.max_nl, noncoll, num_nonloc_ions)
    // init.c or get_ddd get_qqq, we need to check the order
    int proj_index = 0;
    for (int ion = 0; ion < num_nonloc_ions; ion++)
    {

        /*Actual index of the ion under consideration*/
        proj_index = ion * ct.max_nl * ct.max_nl * ct.noncoll_factor * ct.noncoll_factor;
        int gion = nonloc_ions_list[ion];
        SPECIES &AtomType = Species[Atoms[gion].species];

        int nh = AtomType.nh;

        dnmI = Atoms[gion].dnmI;
        qqq = Atoms[gion].qqq;

        for (int i = 0; i < nh; i++)
        {
            int inh = i * nh;
            for (int j = 0; j < nh; j++)
            {

                if(ct.is_ddd_non_diagonal) {
                    if(ct.noncoll)
                    {
                        int it0 = i;
                        int jt0 = j;
                        int it1 = i + ct.max_nl;
                        int jt1 = j + ct.max_nl;
                        M_dnm_C[proj_index + it0 * ct.max_nl * 2 + jt0] = Atoms[gion].dnmI_so[inh+j + 0 * nh *nh];
                        M_dnm_C[proj_index + it0 * ct.max_nl * 2 + jt1] = Atoms[gion].dnmI_so[inh+j + 1 * nh *nh];
                        M_dnm_C[proj_index + it1 * ct.max_nl * 2 + jt0] = Atoms[gion].dnmI_so[inh+j + 2 * nh *nh];
                        M_dnm_C[proj_index + it1 * ct.max_nl * 2 + jt1] = Atoms[gion].dnmI_so[inh+j + 3 * nh *nh];
                        M_qqq_C[proj_index + it0 * ct.max_nl * 2 + jt0] = Atoms[gion].qqq_so[inh+j + 0 * nh *nh];
                        M_qqq_C[proj_index + it0 * ct.max_nl * 2 + jt1] = Atoms[gion].qqq_so[inh+j + 1 * nh *nh];
                        M_qqq_C[proj_index + it1 * ct.max_nl * 2 + jt0] = Atoms[gion].qqq_so[inh+j + 2 * nh *nh];
                        M_qqq_C[proj_index + it1 * ct.max_nl * 2 + jt1] = Atoms[gion].qqq_so[inh+j + 3 * nh *nh];
                    }
                    else
                    {
                        int idx = proj_index + i * ct.max_nl + j;
                        M_dnm[idx] = (KpointType)dnmI[inh+j];
                        M_qqq[idx] = (KpointType)qqq[inh+j];
                    }
                }
                else {
                    // Diagonal for norm conserving so just save those.
                    if(i == j ) {
                        M_dnm[ion*ct.max_nl + j ] = (KpointType)dnmI[inh+j ];
                        M_qqq[ion*ct.max_nl + j ] = (KpointType)qqq[inh+j ];
                    }
                }

            }
        }
    }

    // calculating nvwork_ion (ct.max_nl *noncoll, num_states, ion) =  
    //               M_dnm(ct.max_nl * noncoll, ct.max_nl * noncoll, ion) 
    //               * sint_compack(ct.max_nl*noncoll, num_states, ion)
    //
    if(ct.is_ddd_non_diagonal)
    {
        int dim_a = ct.max_nl * ct.noncoll_factor;
        size_t strideA = (size_t)dim_a * (size_t)dim_a;
        size_t strideB = (size_t)dim_a * (size_t)num_states;
        size_t strideC = (size_t)dim_a * (size_t)num_states;;

        RmgGemmStridedBatched(transa, transa, dim_a, num_states, dim_a, ONE_t, 
                M_dnm, dim_a, strideA, sint_compack, dim_a, strideB, ZERO_t,
                nwork_ion, dim_a, strideC, num_nonloc_ions); 

        //      nvwork_ion: (ct.max_nl, noncoll, num_states, ion)
        // rotate it to nwork (ct.max_nl, ion, noncoll, num_states)`

        for(int st = 0; st < num_states * ct.noncoll_factor; st++)
        {
            for(int ion = 0; ion < num_nonloc_ions; ion++) 
            {
                for(int ih = 0; ih < ct.max_nl; ih++)
                {
                    nwork[st*num_nonloc_ions * ct.max_nl + ion * ct.max_nl + ih]=
                        nwork_ion[ion * strideC + st * ct.max_nl + ih];

                }
            }

        }
    }
    else
    {
        // Optimize for GPU's!
        //for spin-orbit couplling, it always has non_diagonal part.

        for(int jj = 0;jj < num_states;jj++) {
            for (int ion = 0; ion < num_nonloc_ions; ion++) {
                for(int ih = 0;ih < ct.max_nl;ih++) {
                    nwork[jj*num_tot_proj + ion * ct.max_nl + ih] = 
                        M_dnm[ion * ct.max_nl + ih] * 
                        sint_compack[ion * num_states * ct.max_nl + ct.max_nl * jj + ih];
                }
            }
        }
    }


    int tot_states = num_states * ct.noncoll_factor;

    //nwork: num_tot_proj * (ct.noncoll_factor * num_states)

    RmgGemm (transa, transa, P0_BASIS, tot_states, num_tot_proj,
            ONE_t, weight,  P0_BASIS, nwork, num_tot_proj,
            ZERO_t,  nv, P0_BASIS);

    if(! (ct.norm_conserving_pp && ct.is_gamma) ) 
        memcpy(ns, psi, stop*sizeof(KpointType));
    if(!ct.norm_conserving_pp) {

        int dim_a = ct.max_nl * ct.noncoll_factor;
        int strideA = dim_a * dim_a;
        int strideB = dim_a * num_states;
        int strideC = dim_a * num_states;;

        RmgGemmStridedBatched(transa, transa, dim_a, num_states, dim_a, ONE_t, 
                M_qqq, dim_a, strideA, sint_compack, dim_a, strideB, ZERO_t,
                nwork_ion, dim_a, strideC, num_nonloc_ions); 

        //      nvwork_ion: (ct.max_nl, noncoll, num_states, ion)
        // rotate it to nwork (ct.max_nl, ion, noncoll, num_states)`

        for(int st = 0; st < num_states * ct.noncoll_factor; st++)
        {
            for(int ion = 0; ion < num_nonloc_ions; ion++) 
            {
                for(int ih = 0; ih < ct.max_nl; ih++)
                {
                    nwork[st*num_nonloc_ions * ct.max_nl + ion * ct.max_nl + ih]=
                        nwork_ion[ion * strideC + st * ct.max_nl + ih];

                }
            }

        }


        RmgGemm (transa, transa, P0_BASIS, tot_states, num_tot_proj, 
                ONE_t, weight,  P0_BASIS, nwork, num_tot_proj,
                ONE_t,  ns, P0_BASIS);

    }



    if(ct.xc_is_hybrid && Functional::is_exx_active())
    {
        for(size_t i = 0; i < stop; i++) nv[i] += ct.exx_fraction * kpoint->vexx[(size_t)first_state*(size_t)P0_BASIS + i];
        //AppExx(kpoint, psi, num_states, kpoint->vexx, &nv[(size_t)first_state*(size_t)P0_BASIS]);

    }

    RmgFreeHost(M_qqq);
    RmgFreeHost(M_dnm);
    RmgFreeHost(nwork);
    RmgFreeHost(nwork_ion);
    RmgFreeHost(sint_compack);


    // Add in ldaU contributions to nv
    if(ct.ldaU_mode == LDA_PLUS_U_SIMPLE)
    {
        kpoint->ldaU->app_vhubbard(nv, kpoint->orbitalsint_local, first_state, num_states);
    }

}


template void AppS<double>(Kpoint<double> *, double *, double *, double *, int, int);
template void AppS<std::complex<double> >(Kpoint<std::complex<double>> *, std::complex<double> *, 
        std::complex<double> *, std::complex<double> *, int, int);
    template <typename KpointType>
void AppS(Kpoint<KpointType> *kpoint, KpointType *sintR,
        KpointType *psi, KpointType *ns,
        int first_state, int num_states)
{

    // Sanity check
    if(num_states > ct.non_local_block_size)
        throw RmgFatalException() << "AppS called with num_states > non_local_block_size in " << __FILE__ << " at line " << __LINE__ << "\n";

    KpointType *weight = kpoint->nl_weight;
#if HIP_ENABLED || CUDA_ENABLED
    weight = kpoint->nl_weight_gpu;
#endif

    //   sintR:  <beta | psi_up, psi_down>, dimension is numProj * 2 * num_states in noncollinear case.
    //   ns : |psi_up, psu_down > + |beta_n> Dnm <beta_m|psi_up, psi_down>, dimension 2 * pbasis

    int P0_BASIS = kpoint->pbasis;
    int num_nonloc_ions = kpoint->BetaProjector->get_num_nonloc_ions();
    int *nonloc_ions_list = kpoint->BetaProjector->get_nonloc_ions_list();
    int num_tot_proj = kpoint->BetaProjector->get_num_tot_proj();

    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);

    char *transa = "n";

    double *qqq;
    KpointType *psintR;
    size_t stop = (size_t)num_states * (size_t)P0_BASIS * (size_t) ct.noncoll_factor;

    if(num_tot_proj == 0)
    {
        for(size_t idx = 0;idx < stop;idx++) ns[idx] = psi[idx];
        return;
    }


    size_t alloc = (size_t)num_tot_proj * (size_t)num_states * ct.noncoll_factor;
    size_t M_cols = 1;
    if(ct.is_ddd_non_diagonal) M_cols = (size_t)num_tot_proj * ct.noncoll_factor;
    size_t alloc1 = (size_t)num_tot_proj * (size_t)M_cols * ct.noncoll_factor;

    KpointType *sint_compack = (KpointType *)RmgMallocHost(sizeof(KpointType) * alloc);
    KpointType *nwork = (KpointType *)RmgMallocHost(sizeof(KpointType) * alloc);
    KpointType *M_qqq = (KpointType *)RmgMallocHost(sizeof(KpointType) * alloc1);
    for(size_t i = 0;i < alloc;i++) sint_compack[i] = 0.0;
    std::complex<double> *M_qqq_C = (std::complex<double> *) M_qqq;

    for(int istate = 0; istate < num_states * ct.noncoll_factor; istate++)
    {
        int sindex = (istate + first_state * ct.noncoll_factor) * num_nonloc_ions * ct.max_nl;
        for (int ion = 0; ion < num_nonloc_ions; ion++)
        {
            int proj_index = ion * ct.max_nl;
            psintR = &sintR[proj_index + sindex];
            //psintI = &sintI[ion * num_states * ct.max_nl + sindex];
            /*Actual index of the ion under consideration*/
            int gion = nonloc_ions_list[ion];
            SPECIES &AtomType = Species[Atoms[gion].species];

            int nh = AtomType.nh;
            for (int i = 0; i < nh; i++)
            {
                sint_compack[istate * num_tot_proj + proj_index + i] = psintR[i];
            }
        }
    }

    for (size_t i = 0; i < alloc1; i++)
    {
        M_qqq[i] = ZERO_t;
    }


    // set up M_qqq this can be done outside in the
    // init.c or get_ddd get_qqq, we need to check the order
    int proj_index = 0;
    for (int ion = 0; ion < num_nonloc_ions; ion++)
    {

        /*Actual index of the ion under consideration*/
        proj_index = ion * ct.max_nl;
        int gion = nonloc_ions_list[ion];
        SPECIES &AtomType = Species[Atoms[gion].species];

        int nh = AtomType.nh;

        qqq = Atoms[gion].qqq;

        for (int i = 0; i < nh; i++)
        {
            int inh = i * nh;
            for (int j = 0; j < nh; j++)
            {

                if(ct.is_ddd_non_diagonal) {
                    if(ct.noncoll)
                    {
                        int it0 = proj_index + i;
                        int jt0 = proj_index + j;
                        int it1 = proj_index + i + num_tot_proj;
                        int jt1 = proj_index + j + num_tot_proj;
                        M_qqq_C[it0 * num_tot_proj * 2 + jt0] = Atoms[gion].qqq_so[inh+j + 0 * nh *nh];
                        M_qqq_C[it0 * num_tot_proj * 2 + jt1] = Atoms[gion].qqq_so[inh+j + 1 * nh *nh];
                        M_qqq_C[it1 * num_tot_proj * 2 + jt0] = Atoms[gion].qqq_so[inh+j + 2 * nh *nh];
                        M_qqq_C[it1 * num_tot_proj * 2 + jt1] = Atoms[gion].qqq_so[inh+j + 3 * nh *nh];
                    }
                    else
                    {
                        int idx = (proj_index + i) * num_tot_proj + proj_index + j;
                        M_qqq[idx] = (KpointType)qqq[inh+j];
                    }
                }
                else {
                    // Diagonal for norm conserving so just save those.
                    if((proj_index + i) == (proj_index + j)) {
                        M_qqq[proj_index + j ] = (KpointType)qqq[inh+j ];
                    }
                }

            }
        }
    }



    int dim_dnm = num_tot_proj * ct.noncoll_factor;
    int tot_states = num_states * ct.noncoll_factor;
    if(!ct.norm_conserving_pp) {

        //sint_compack: dim_dnm * num_states == num_tot_proj * ct.noncoll_factor * num_states
        //nwork: dim_dnm * num_states == num_tot_proj * ct.noncoll_factor * num_states
        //  in the first RmgGemm, nwork is a matrix of (dim_dnm) * num_states 
        //  in the second RmgGemm, nwork is a matrix of num_tot_proj * (tot_states) 

        // leading dimension is num_tot_proj * 2 for noncollinear
        memcpy(ns, psi, stop*sizeof(KpointType));

        RmgGemm (transa, transa, dim_dnm, num_states, dim_dnm, 
                ONE_t, M_qqq,  dim_dnm, sint_compack, dim_dnm,
                ZERO_t,  nwork, dim_dnm);

        RmgGemm (transa, transa, P0_BASIS, tot_states, num_tot_proj, 
                ONE_t, weight,  P0_BASIS, nwork, num_tot_proj,
                ONE_t,  ns, P0_BASIS);

    }
    else 
    {

#if CUDA_ENABLED
        // For norm conserving and gamma ns=psi so other parts of code were updated to not require this
        gpuMemcpy(ns, psi, stop*sizeof(KpointType), gpuMemcpyDefault);
#else
        memcpy(ns, psi, stop*sizeof(KpointType));
#endif

    }

    RmgFreeHost(M_qqq);
    RmgFreeHost(nwork);
    RmgFreeHost(sint_compack);

}

template <typename T> void AppExx(Kpoint<double> *, double *, int, double *, double);
template <typename T> void AppExx(Kpoint<std::complex<double>> *, std::complex<double> *, int, std::complex<double> *, std::complex<double> *);

template <typename T> void AppExx(Kpoint<T> *kptr, T *psi, int N, T *vexx, T *nv)
{
    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a = trans_t;
    int pbasis = kptr->pbasis;
    int factor=1;
    if(typeid(T) == typeid(std::complex<double>)) factor = 2;
    if(typeid(T) == typeid(std::complex<double>)) trans_a = trans_c;
    double vel = kptr->L->get_omega() / ((double)(kptr->G->get_NX_GRID(1) * kptr->G->get_NY_GRID(1) * kptr->G->get_NZ_GRID(1)));
    T alpha(1.0);
    T alphavel(vel);
    T beta(0.0);
    T exx_fraction(ct.exx_fraction);

    // Compute the overlap matrix
    T *overlaps = new T[kptr->nstates * N];

    RmgGemm(trans_a, trans_n, N, kptr->nstates, pbasis, alphavel, kptr->prev_orbitals, pbasis,
            psi, pbasis, beta, overlaps, N);

    BlockAllreduce((double *)overlaps, (size_t)(kptr->nstates)*(size_t)N * (size_t)factor, kptr->grid_comm);

    // Update nv
    RmgGemm(trans_n, trans_n, pbasis, N, kptr->nstates, exx_fraction, vexx, pbasis,
            overlaps, kptr->nstates, alpha, nv, pbasis);

    delete [] overlaps;
}

