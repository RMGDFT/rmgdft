/*
 *
 * Copyright 2025 The RMG Project Developers. See the COPYRIGHT file 
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

#include <complex>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "Gpufuncs.h"
#include "Subdiag.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"
#include "Solvers.h"
#include "Functional.h"
#include "RmgMatrix.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"


extern Kpoint<double> **Kptr_g;

void compute_vxc(spinobj<double> &rho, fgobj<double> &rhocore, double &XC, double &vtxc, spinobj<double> &v_xc, int nspin)
{
    compute_vxc(rho.data(), rhocore.data(), XC, vtxc, v_xc.data(), nspin);
}

void compute_vxc(double *rho, double *rhocore, double &XC, double &vtxc, double *v_xc, int nspin)
{
    static Functional *F = NULL;
    if(F == NULL) F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, nspin);

    if(F->dft_is_meta_rmg())
    {
        wfobj<double> kdetau_c;
        fgobj<double> kdetau_f;
        kdetau_c.set(0.0);
        if(ct.scf_steps >= 0)
        {
            for(int ik = 0; ik < ct.num_kpts_pe; ik++) Kptr_g[ik]->KineticEnergyDensity(kdetau_c.data());
//            FftInterpolation(*Rmg_G, kdetau_c.data(), kdetau_f.data(), 2, false);
            int ratio = Rmg_G->default_FG_RATIO;
            Prolong P(2, ct.prolong_order, 0.0, *Rmg_T,  Rmg_L, *Rmg_G);
            int dimx = kdetau_f.dimx;
            int dimy = kdetau_f.dimy;
            int dimz = kdetau_f.dimz;
            int half_dimx = kdetau_c.dimx;
            int half_dimy = kdetau_c.dimy;
            int half_dimz = kdetau_c.dimz;
            P.prolong(kdetau_f.data(), kdetau_c.data(), dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);

        }
//  Need to update GetNewRho if we want to compute this on fine grid but not clear if it's necessary
//        else
//        {
//            for(int ix=0;ix < this->pbasis;ix++) kdetau_f[ix] = this->ke_density[ix];
//        }
        for(int ix = 0;ix < nspin*kdetau_f.pbasis;ix++) v_xc[ix] = 0.0;
        F->v_xc_meta(rho, rhocore, XC, vtxc, v_xc, kdetau_f.data(), nspin);
        return;
    }

    F->v_xc(rho, rhocore, XC, vtxc, v_xc, nspin );
//    delete F;
}
