/*
 *
 * Copyright 2026 The RMG Project Developers. See the COPYRIGHT file 
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
#include <fcntl.h>
#ifdef USE_NUMA
    #include <numa.h>
#endif
#include <sys/mman.h>
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
#include "Solvers.h"
#include "Subdiag.h"
#include "Functional.h"
#include "GpuAlloc.h"
#include "rmg_sum_all.h"

#include "RmgException.h"
#include "Functional.h"
#include "FiniteDiff.h"
#include "rmgthreads.h"


/* this function is used to approximate the spectral radius of the preconditioned
   Hamiltonian. It is a rough approximation since it doesn't account for the
   non-local potential operators.
*/

template double spectral_radius<double>(fgobj<double> &, spinobj<double> &, Kpoint<double>*);
template double spectral_radius<std::complex<double> >(fgobj<double> &, spinobj<double> &, Kpoint<std::complex <double> >*);

template <typename OrbitalType> double spectral_radius (fgobj<double> &vtot, spinobj<double> &vxc, Kpoint<OrbitalType> *kptr)
{
    wfobj<double> vtot_psi;
    double hxgrid = Rmg_G->get_hxgrid(1);
    FiniteDiff FD(&Rmg_L, ct.alt_laplacian);
    double diag = FD.fd_coeff0(ct.kohn_sham_fd_order, hxgrid);

    /*dnmI has to be stup before calling subdiag */
    double *vxc_psi = NULL;

    // Transfer vtot from the fine grid to the wavefunction grid for Subdiag
    GetVtotPsi (vtot_psi.data(), vtot.data(), Rmg_G->default_FG_RATIO);
    if(ct.noncoll)
    {
        vxc_psi = new double[vtot.pbasis];
        for(int is = 0; is < 4; is++)
            GetVtotPsi (&vxc_psi[is*vtot_psi.pbasis], &vxc[is*vtot.pbasis], Rmg_G->default_FG_RATIO);
    }
    State<OrbitalType> *sp = &kptr->Kstates[0];
    wfobj<OrbitalType> psi0, u, Hu;
    for(int i=0;i < psi0.pbasis;i++) psi0[i] = sp->psi[i];    
    long int idum = (long int)pct.gridpe;
    rand0 (&idum);
    for (int i=0;i < u.pbasis;i++) sp->psi[i] = rand0 (&idum) - 0.5;

    double eig;
    for(int its=0;its < 20;its++)
    {
        kptr->BetaProjector->project(kptr, kptr->newsint_local, 0, 
                0, kptr->nl_weight);
        AppNls(kptr, kptr->newsint_local, sp->psi, kptr->nv, kptr->ns, 0, 1);
        ApplyHamiltonian<OrbitalType,OrbitalType> (kptr, sp, sp->istate, sp->psi, Hu.data(),
                     vtot_psi.data(), vxc_psi, kptr->nv, false);
        for(int i=0;i < u.pbasis;i++) Hu[i] /= (-0.5*diag + vtot_psi[i]);
        eig = ComputeEig(u.pbasis, sp->psi, Hu.data(), sp->psi);
        for(int i=0;i < u.pbasis;i++) sp->psi[i] = Hu[i];
        double sum = 0.0;
        for(int i=0;i < u.pbasis;i++) sum += std::norm(sp->psi[i])*get_vel();
        sum = rmg::sum_all(sum, pct.grid_comm);
        sum = 1.0 / sqrt(sum);
        for(int i=0;i < u.pbasis;i++) sp->psi[i] *= sum;
    }
    //if(pct.gridpe==0)printf("EEEE  %14.8f\n",eig);
    for(int i=0;i < psi0.pbasis;i++) kptr->Kstates[0].psi[i] = psi0[i];    
    if(vxc_psi) delete [] vxc_psi;
    return eig;
}



