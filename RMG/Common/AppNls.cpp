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
#include "../Headers/prototypes.h"

#include "GlobalSums.h"
#include "blas.h"

template void AppNls<double>(Kpoint<double> *, double *);
template void AppNls<std::complex<double> >(Kpoint<std::complex<double>> *, std::complex<double> *);


template <typename KpointType>
void AppNls(Kpoint<KpointType> *kpoint, KpointType *sintR)
{

    int num_states = kpoint->get_nstates();
    int P0_BASIS = kpoint->pbasis;
    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);
    KpointType *NULLptr = NULL;

    char *transa = "n";

    double *dnmI;
    double *qqq;
    KpointType *psintR;

    KpointType *psi = kpoint->orbital_storage;
    KpointType *nv = kpoint->nv;
    KpointType *ns = kpoint->ns;
    KpointType *Bns = kpoint->Bns;

    if(pct.num_tot_proj == 0)
    {
        for(int i = 0; i < num_states * P0_BASIS; i++) nv[i] = ZERO_t;
        if(!ct.norm_conserving_pp) for(int i = 0; i < num_states * P0_BASIS; i++) Bns[i] = ZERO_t;
        for(int idx = 0;idx < num_states * P0_BASIS;idx++) ns[idx] = psi[idx];

        return;
    }


    int alloc = pct.num_tot_proj * num_states;
#if GPU_ENABLED
    KpointType *sint_compack = (KpointType *)GpuMallocHost(sizeof(KpointType) * alloc);
    KpointType *nwork = (KpointType *)GpuMallocHost(sizeof(KpointType) * alloc);
    KpointType *M_dnm = (KpointType *)GpuMallocHost(sizeof(KpointType) * pct.num_tot_proj * pct.num_tot_proj);
    KpointType *M_qqq = (KpointType *)GpuMallocHost(sizeof(KpointType) * pct.num_tot_proj * pct.num_tot_proj);
    for(int i = 0;i < alloc;i++) sint_compack[i] = 0.0;
#else
    KpointType *sint_compack = new KpointType[alloc]();
    KpointType *nwork = new KpointType[alloc];
    KpointType *M_dnm = new KpointType [pct.num_tot_proj * pct.num_tot_proj];
    KpointType *M_qqq = new KpointType [pct.num_tot_proj * pct.num_tot_proj];
#endif


    for(int istate = 0; istate < num_states; istate++)
    {
        int sindex = istate * ct.max_nl;
        for (int ion = 0; ion < pct.num_nonloc_ions; ion++)
        {
            int proj_index = ion * ct.max_nl;
            psintR = &sintR[ion * num_states * ct.max_nl + sindex];
            //psintI = &sintI[ion * num_states * ct.max_nl + sindex];
            /*Actual index of the ion under consideration*/
            int gion = pct.nonloc_ions_list[ion];
            ION *iptr = &ct.ions[gion];
            SPECIES *sp = &ct.sp[iptr->species];

            int nh = sp->nh;
            for (int i = 0; i < nh; i++)
            {
                sint_compack[istate * pct.num_tot_proj + proj_index + i] = psintR[i];
            }
        }
    }

    for (int i = 0; i < pct.num_tot_proj * pct.num_tot_proj; i++)
    {
        M_dnm[i] = ZERO_t;
        M_qqq[i] = ZERO_t;
        pct.M_dnm[i] = 0.0;
        pct.M_qqq[i] = 0.0;
    }


    // set up pct.M_qqq and pct.M_dnm, this can be done outside in the
    // init.c or get_ddd get_qqq, we need to check the order
    int proj_index = 0;
    for (int ion = 0; ion < pct.num_nonloc_ions; ion++)
    {

        /*Actual index of the ion under consideration*/
        proj_index = ion * ct.max_nl;
        int gion = pct.nonloc_ions_list[ion];
        ION *iptr = &ct.ions[gion];
        SPECIES *sp = &ct.sp[iptr->species];

        int nh = sp->nh;

        dnmI = pct.dnmI[gion];
        qqq = pct.qqq[gion];

        for (int i = 0; i < nh; i++)
        {
            int inh = i * nh;
            for (int j = 0; j < nh; j++)
            {

                int idx = (proj_index + i) * pct.num_tot_proj + proj_index + j;
                M_dnm[idx] = (KpointType)dnmI[inh+j];
                M_qqq[idx] = (KpointType)qqq[inh+j];
                pct.M_dnm[idx] = std::real(dnmI[inh+j]);
                pct.M_qqq[idx] = std::real(qqq[inh+j]);

            }
        }
    }

    RmgGemm (transa, transa, pct.num_tot_proj, num_states, pct.num_tot_proj, 
            ONE_t, M_dnm,  pct.num_tot_proj, sint_compack, pct.num_tot_proj,
            ZERO_t,  nwork, pct.num_tot_proj, NULLptr, NULLptr, NULLptr, false, false, false, true);


    RmgGemm (transa, transa, P0_BASIS, num_states, pct.num_tot_proj, 
            ONE_t, kpoint->nl_Bweight,  P0_BASIS, nwork, pct.num_tot_proj,
            ZERO_t,  nv, P0_BASIS, NULLptr, NULLptr, NULLptr, false, false, false, true);


    for(int idx = 0;idx < num_states * P0_BASIS;idx++)
        ns[idx] = psi[idx];

    if(!ct.norm_conserving_pp) {

        RmgGemm (transa, transa, pct.num_tot_proj, num_states, pct.num_tot_proj, 
                ONE_t, M_qqq,  pct.num_tot_proj, sint_compack, pct.num_tot_proj,
                ZERO_t,  nwork, pct.num_tot_proj, NULLptr, NULLptr, NULLptr, false, false, false, true);

        RmgGemm (transa, transa, P0_BASIS, num_states, pct.num_tot_proj, 
                ONE_t, kpoint->nl_weight,  P0_BASIS, nwork, pct.num_tot_proj,
                ONE_t,  ns, P0_BASIS, NULLptr, NULLptr, NULLptr, false, false, false, true);

        RmgGemm (transa, transa, P0_BASIS, num_states, pct.num_tot_proj, 
                ONE_t, kpoint->nl_Bweight,  P0_BASIS, nwork, pct.num_tot_proj,
                ZERO_t,  Bns, P0_BASIS, NULLptr, NULLptr, NULLptr, false, false, false, true);
    }

#if GPU_ENABLED
    GpuFreeHost(M_qqq);
    GpuFreeHost(M_dnm);
    GpuFreeHost(nwork);
    GpuFreeHost(sint_compack);
#else
    delete [] nwork;
    delete [] sint_compack;
    delete [] M_qqq;
    delete [] M_dnm;
#endif

}

