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
#include "Solvers.h"

template void BandStructure(Kpoint<double> **, double *vh, double *vxc, double *vnuc, std::vector<bool> exclude_bands);
template void BandStructure(Kpoint<std::complex<double> > **, double *vh, double *vxc, double *vnuc, std::vector<bool> exclude_bands);

template <typename KpointType>
void BandStructure(Kpoint<KpointType> ** Kptr, double *vh, double *vxc, double *vnuc, std::vector<bool> exclude_bands)
{
    
    double *vtot, *vtot_psi, *vxc_psi=NULL;
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

    get_ddd (vtot, vxc, true);
    GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);

    if(ct.noncoll)
    {   
        vxc_psi = new double[ct.nspin * P0_BASIS];
        for(int is = 0; is < ct.nspin; is++)
            GetVtotPsi (&vxc_psi[is*P0_BASIS], &vxc[is*FP0_BASIS], Rmg_G->default_FG_RATIO);
    }

    if(ct.alloc_states < 2 * ct.num_states) 
    {
        printf("\n band strucuture need 2*ct.num_states allocated\n");
        exit(0);
    }

    // Broken by Lcao changes
    //    LcaoGetPsi(&Kptr[0]->Kstates[ct.num_states]);

    int projector_type = DELOCALIZED;
    if(ct.localize_projectors) projector_type = LOCALIZED;

    // Number of projectors required is computed when the Projector is created.
    // Beta function weights are created in the calls to get_nlop.
    Kptr[0]->get_nlop(projector_type);
    for(int kpt=1; kpt < ct.num_kpts_pe; kpt++)
    {
        if(Kptr[kpt]->BetaProjector) delete Kptr[kpt]->BetaProjector;

        Kptr[kpt]->BetaProjector = new Projector<KpointType>(projector_type, ct.max_nl, BETA_PROJECTOR);
        Kptr[kpt]->nl_weight_size = (size_t)Kptr[kpt]->BetaProjector->get_num_tot_proj() * (size_t)Kptr[kpt]->pbasis + 128;
        Kptr[kpt]->nl_weight = Kptr[0]->nl_weight;
        Kptr[kpt]->newsint_local = Kptr[0]->newsint_local;

    } // end loop over kpts



    // Loop over k-points
    double *res = new double[ct.num_states];
    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++) {

        if(ct.localize_projectors)
        {
            Kptr[kpt]->GetLocalizedWeight ();
        }
        else
        {
            Kptr[kpt]->GetDelocalizedWeight ();
        }

        if((ct.ldaU_mode != LDA_PLUS_U_NONE) && (ct.num_ldaU_ions > 0))
        {
            Kptr[kpt]->GetDelocalizedOrbital ();
            Kptr[kpt]->get_ldaUop(ct.atomic_orbital_type);
        }

        Kptr[kpt]->nstates = ct.num_states;

        if (Verify ("kohn_sham_solver","multigrid", Kptr[0]->ControlMap) ) {

            for (ct.scf_steps = 0, CONVERGED = false;
                    ct.scf_steps < ct.max_scf_steps && !CONVERGED; ct.scf_steps++)
            {
                RmgTimer *RT = new RmgTimer("MgridSub in band");

                //Kptr[kpt]->LcaoGetPsi(ct.num_states, ct.num_states);
                Kptr[kpt]->MgridSubspace(vtot_psi, vxc_psi);

                delete RT;

                for(int istate = 0; istate < ct.num_states; istate++)
                    res[istate] = Kptr[kpt]->Kstates[istate].res;
                GlobalSums(res, ct.num_states, pct.grid_comm);


                max_res = res[0];
                for(int istate = 0; istate < ct.num_states; istate++)
                    max_res = std::max(max_res, res[istate]);
                //        for(int istate = 0; istate < Kptr[kpt]->nstates; istate++)
                //         rmg_printf("\n kpt = %d scf =%d state=%d res = %e", kpt, ct.scf_steps, istate, Kptr[kpt]->Kstates[istate].res);
                if(pct.gridpe == 0) printf("\n kpt= %d  scf = %d  max_res = %e", kpt+pct.kstart, ct.scf_steps, max_res);

                if (max_res <1.0e-3)
                {
                    rmg_printf("\n BAND STRUCTURE: Converged with max_res %10.5e", max_res);
                    CONVERGED = true;
                }

            } // end loop scf
        }
        else if(Verify ("kohn_sham_solver","davidson", Kptr[0]->ControlMap)) {
            int notconv;
            ct.scf_steps = 2;
            //   ct.scf_accuracy = 1.0e-10;
            Kptr[kpt]->Davidson(vtot_psi, vxc_psi, notconv);
        }

        //for(int istate = 0; istate < Kptr[kpt]->nstates; istate++)
        //rmg_printf("\n BAND STRUCTURE: state %d res %10.5e ", istate, Kptr[kpt]->Kstates[istate].res);

        int kpt_global = kpt + pct.kstart;
        if(ct.wannier90) {
            int factor = 1;
            if(ct.nspin == 2) factor = 2;
            ct.kp[kpt+pct.kstart].eigs.resize(ct.num_states * factor);
            for(int st = 0; st < ct.num_states; st++)
                ct.kp[kpt+pct.kstart].eigs[st] = Kptr[kpt]->Kstates[st].eig[0];

            if(ct.nspin == 2)
            {
                for(int st = 0; st < ct.num_states; st++)
                    ct.kp[kpt+pct.kstart].eigs[st + ct.num_states] = Kptr[kpt]->Kstates[st].eig[1];
            }

            Write_Wfs_forWannier(kpt_global, Kptr[kpt], exclude_bands, "WfsForWannier90/wfs");
        }
        if(ct.rmg2bgw)
        {
            WriteBGW_Wfng(kpt_global, Kptr[kpt]);
            WriteBGW_VxcEig(kpt_global,vxc, Kptr[kpt]);
        }


    } // end loop over kpoints


    delete [] res;
    delete [] vtot;
    delete [] vtot_psi;
    if(ct.noncoll) delete [] vxc_psi;
}

/******/

