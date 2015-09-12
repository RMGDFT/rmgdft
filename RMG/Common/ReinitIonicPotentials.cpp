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

#include "make_conf.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Kpoint.h"
#include "ErrorFuncs.h"
#include "transition.h"

template void ReinitIonicPotentials<double>(Kpoint<double> **, double *, double *, double *);
template void ReinitIonicPotentials<std::complex<double> >(Kpoint<std::complex<double>> **, double *, double *, double *);

template <typename KpointType>
void ReinitIonicPotentials (Kpoint<KpointType> **Kptr, double * vnuc, double * rhocore, double * rhoc)
{

    int pbasis = Kptr[0]->pbasis;

    /* Update items that change when the ionic coordinates change */
    init_nuc (vnuc, rhoc, rhocore);
    get_QI ();
    GetNlop(Kptr);

    // Number of total projectors required is computed in GetNlop so we allocate per
    // k-point storage for the weights here.
    for(int kpt=0;kpt < ct.num_kpts;kpt++) {

        if(ct.is_gamma) {

            // Identical for gamma point
            Kptr[kpt]->nl_weight = (KpointType *)pct.weight;
            Kptr[kpt]->nl_Bweight = (KpointType *)pct.Bweight;

        }
        else {

#if GPU_ENABLED
            if(Kptr[kpt]->nl_weight != NULL) cudaFreeHost(Kptr[kpt]->nl_weight);
            if(Kptr[kpt]->nl_Bweight != NULL) cudaFreeHost(Kptr[kpt]->nl_Bweight);

            cudaError_t cuerr;
            // Allocate new storage
            if(pct.num_tot_proj) {

                cuerr = cudaMallocHost((void **)&Kptr[kpt]->nl_weight , pbasis * pct.num_tot_proj * sizeof(KpointType) );
                if(cuerr != cudaSuccess)
                    RmgCudaError(__FILE__, __LINE__, cuerr, "GPU memory allocation error");

             cuerr = cudaMallocHost((void **)&Kptr[kpt]->nl_Bweight , pbasis * pct.num_tot_proj * sizeof(KpointType) );
                if(cuerr != cudaSuccess)
                    RmgCudaError(__FILE__, __LINE__, cuerr, "GPU memory allocation error");

            }
#else
            // Release old memory storage for weights
            if(Kptr[kpt]->nl_weight != NULL) delete [] Kptr[kpt]->nl_weight;
            if(Kptr[kpt]->nl_Bweight != NULL) delete [] Kptr[kpt]->nl_Bweight;

            // Allocate new storage
            if(pct.num_tot_proj) {
                Kptr[kpt]->nl_weight = new KpointType[pct.num_tot_proj * pbasis]();
                Kptr[kpt]->nl_Bweight = new KpointType[pct.num_tot_proj * pbasis]();
            }
#endif
        }

    } // end loop over kpts


    /*Other things that need to be recalculated when ionic positions change */
    GetWeight (Kptr);
    get_qqq ();


#if 0
    if (!verify ("calculation_mode", "Band Structure Only"))
    {
        betaxpsi (states);
        mix_betaxpsi(0);
    }
#endif

}                               /* end reinit_ionic_pp */

