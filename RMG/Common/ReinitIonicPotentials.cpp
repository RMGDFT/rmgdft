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


#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Kpoint.h"
#include "Atomic.h"
#include "ErrorFuncs.h"
#include "GpuAlloc.h"
#include "transition.h"

template void ReinitIonicPotentials<double>(Kpoint<double> **, double *, double *, double *);
template void ReinitIonicPotentials<std::complex<double> >(Kpoint<std::complex<double>> **, double *, double *, double *);

template <typename KpointType>
void ReinitIonicPotentials (Kpoint<KpointType> **Kptr, double * vnuc, double * rhocore, double * rhoc)
{
    RmgTimer *RT1;
    RmgTimer RT0("3-ReinitIonicPotentials");
    int pbasis = Kptr[0]->pbasis;
    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);


    /* Update items that change when the ionic coordinates change */
    RT1= new RmgTimer("3-ReinitIonicPotentials: init_nuc");
    //init_nuc (vnuc, rhoc, rhocore);
    double *dum_array = NULL;
    if(ct.localize_localpp) 
    {
        InitLocalObject (vnuc, dum_array, ATOMIC_LOCAL_PP, false);
        InitLocalObject (rhoc, dum_array, ATOMIC_RHOCOMP, false);
        InitLocalObject (rhocore, dum_array, ATOMIC_RHOCORE, false);
    }
    else
    {
        InitDelocalizedObject(vnuc, dum_array, ATOMIC_LOCAL_PP, false);
        InitDelocalizedObject (rhocore, dum_array, ATOMIC_RHOCORE, false);
        // For delocalized rhoc is zero
        for(int ix=0;ix < FP0_BASIS;ix++) rhoc[ix] = 0.0;
    }

    //    InitLocalBackward(vnuc, rhoc, rhocore);
    delete RT1;
    RT1= new RmgTimer("3-ReinitIonicPotentials: get_QI");
    GetQI ();
    delete RT1;
    RT1= new RmgTimer("3-ReinitIonicPotentials: GetNlop");
    GetNlop(Kptr);
    delete RT1;

    // Number of total projectors required is computed in GetNlop so we allocate per
    // k-point storage for the weights here.
    size_t alloc = (size_t)pbasis * (size_t)pct.num_tot_proj;
    for(int kpt=0; kpt < ct.num_kpts_pe; kpt++)
    {

        if(ct.is_gamma) {

            // Identical for gamma point
            Kptr[kpt]->nl_weight = (KpointType *)pct.weight;
            Kptr[kpt]->nl_Bweight = (KpointType *)pct.Bweight;
        }
        else {

#if GPU_ENABLED
            if(Kptr[kpt]->nl_weight != NULL) GpuFreeManaged(Kptr[kpt]->nl_weight);
            if((Kptr[kpt]->nl_Bweight != NULL) && ct.need_Bweight) GpuFreeManaged(Kptr[kpt]->nl_Bweight);

            // Allocate new storage
            if(pct.num_tot_proj) 
            {

                Kptr[kpt]->nl_weight = (KpointType *)GpuMallocManaged(alloc * sizeof(KpointType));
                if(ct.need_Bweight) 
                {
                    Kptr[kpt]->nl_Bweight = (KpointType *)GpuMallocManaged(alloc * sizeof(KpointType));
                }
                else 
                {
                    Kptr[kpt]->nl_Bweight = Kptr[kpt]->nl_weight;
                }

            }
#else
            // Release old memory storage for weights
            if(Kptr[kpt]->nl_weight != NULL) delete [] Kptr[kpt]->nl_weight;
            if((Kptr[kpt]->nl_Bweight != NULL) && ct.need_Bweight) delete [] Kptr[kpt]->nl_Bweight;

            // Allocate new storage
            if(pct.num_tot_proj) 
            {
                Kptr[kpt]->nl_weight = new KpointType[alloc]();
                if(ct.need_Bweight) 
                {
                    Kptr[kpt]->nl_Bweight = new KpointType[alloc]();
                }
                else
                {
                    Kptr[kpt]->nl_Bweight = Kptr[kpt]->nl_weight;
                }

            }
#endif
        }

    } // end loop over kpts


    /*Other things that need to be recalculated when ionic positions change */
    RT1= new RmgTimer("3-ReinitIonicPotentials: GetWeight");
    if(ct.localize_projectors) {
        GetWeightLocal (Kptr);
    }
    else {
        GetWeight (Kptr);
    }
    delete RT1;

    RT1= new RmgTimer("3-ReinitIonicPotentials: get_qqq");
    get_qqq ();
    delete RT1;


}                               /* end ReinitIonicPotentials */

