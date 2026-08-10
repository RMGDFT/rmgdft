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
#include "rmg_reduce.h"
#include "Kpoint.h"
#include "rmg_gemm.h"
#include "Gpufuncs.h"
#include "Subdiag.h"
#include "GpuAlloc.h"

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
    Functional *F = NULL;
    //if(F == NULL) F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, nspin);
    F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, nspin);

    if(F->dft_is_meta_rmg())
    {
        wfobj<double> kdetau_c;
        fgobj<double> kdetau_f;
        kdetau_c.set(0.0);
        for(int ix=0;ix < kdetau_f.pbasis;ix++) kdetau_f[ix] = F->ke_density[ix];
        for(int ix = 0;ix < nspin*kdetau_f.pbasis;ix++) v_xc[ix] = 0.0;
        F->v_xc_meta(rho, rhocore, XC, vtxc, v_xc, kdetau_f.data(), nspin);
    }
    else
    {
        F->v_xc(rho, rhocore, XC, vtxc, v_xc, nspin );
    }
    delete F;
}
