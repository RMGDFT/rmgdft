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
#include <stdlib.h>
#include <stdio.h>
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "BaseThread.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "RmgException.h"
#include "GlobalSums.h"
#include "rmgthreads.h"
#include "vhartree.h"
#include "packfuncs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "../Headers/prototypes.h"
#include "Solvers.h"

template void BandStructure(Kpoint<double> **, double *vh, double *vxc, double *vnuc);
template void BandStructure(Kpoint<std::complex<double> > **, double *vh, double *vxc, double *vnuc);

template <typename KpointType>
void BandStructure(Kpoint<KpointType> ** Kptr, double *vh, double *vxc, double *vnuc)
{
    
    double *vtot, *vtot_psi;
    double max_res;
    int FP0_BASIS, P0_BASIS;
    
    int idx;

    bool CONVERGED;

    
    P0_BASIS =  Rmg_G->get_P0_BASIS(1);

    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    vtot = new  double[FP0_BASIS];
    vtot_psi = new double[P0_BASIS];

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];

    get_ddd (vtot);
    GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);

    if(ct.alloc_states < 2 * ct.num_states) 
    {
        printf("\n band strucuture need 2*ct.num_states allocated\n");
        exit(0);
    }

    LcaoGetPsi(&Kptr[0]->Kstates[ct.num_states]);


    // Loop over k-points
    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++) {

        Kptr[kpt]->nstates = 2 * ct.num_states;

        for (ct.scf_steps = 0, CONVERGED = false;
                ct.scf_steps < ct.max_scf_steps && !CONVERGED; ct.scf_steps++)
        {
            RmgTimer *RT = new RmgTimer("MgridSub in band");
            MgridSubspace(Kptr[kpt], vtot_psi);
            delete RT;


            max_res = Kptr[kpt]->Kstates[0].res;
            for(int istate = 0; istate < ct.num_states; istate++)
                if( max_res < Kptr[kpt]->Kstates[istate].res) 
                    max_res = Kptr[kpt]->Kstates[istate].res;
            //         for(int istate = 0; istate < Kptr[kpt]->nstates; istate++)
             //          rmg_printf("\n kpt = %d scf =%d state=%d res = %e", kpt, ct.scf_steps, istate, Kptr[kpt]->Kstates[istate].res);
            //if(pct.gridpe == 0) printf("\n kpt= %d  scf = %d  max_res = %e", kpt+pct.kstart, ct.scf_steps, max_res);

            //if (max_res <1.0e-3)
           // {
            //    rmg_printf("\n BAND STRUCTURE: Converged with max_res %10.5e", max_res);
             //   CONVERGED = true;
           // }

        } // end loop scf

        //for(int istate = 0; istate < Kptr[kpt]->nstates; istate++)
        //rmg_printf("\n BAND STRUCTURE: state %d res %10.5e ", istate, Kptr[kpt]->Kstates[istate].res);

        int kpt_global = kpt + pct.kstart;
        if(ct.rmg2bgw)
        {
            WriteBGW_Wfng(kpt_global, Kptr[kpt]);
            WriteBGW_VxcEig(kpt_global,vxc, Kptr[kpt]);
        }


    } // end loop over kpoints


    delete [] vtot;
    delete [] vtot_psi;
}

/******/





