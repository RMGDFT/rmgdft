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

    delete RT1;
    RT1= new RmgTimer("3-ReinitIonicPotentials: get_QI");
    GetQI ();
    delete RT1;

    int projector_type = DELOCALIZED;
    if(ct.localize_projectors) projector_type = LOCALIZED;
    RT1= new RmgTimer("3-ReinitIonicPotentials: GetNlop");

    // Number of projectors required is computed when the Projector is created.
    // Beta function weights are created in the calls to get_nlop.
    for(int kpt=0; kpt < ct.num_kpts_pe; kpt++)
    {
        Kptr[kpt]->get_nlop(projector_type);
        if((ct.ldaU_mode != LDA_PLUS_U_NONE) && (ct.num_ldaU_ions > 0))
        {
            Kptr[kpt]->get_ldaUop(ct.atomic_orbital_type);
        }

    } // end loop over kpts
    delete RT1;


    /*Other things that need to be recalculated when ionic positions change */
    RT1= new RmgTimer("3-ReinitIonicPotentials: GetWeight");
    if(ct.localize_projectors) {
        GetLocalizedWeight (Kptr);
    }
    else {
        GetDelocalizedWeight (Kptr);
    }

    if((ct.ldaU_mode != LDA_PLUS_U_NONE) && (ct.num_ldaU_ions > 0))
    {
        GetDelocalizedOrbital (Kptr);
    }

    delete RT1;

    RT1= new RmgTimer("3-ReinitIonicPotentials: get_qqq");
    get_qqq ();
    delete RT1;


}                               /* end ReinitIonicPotentials */

