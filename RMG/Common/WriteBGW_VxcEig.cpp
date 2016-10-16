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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#if !(defined(_WIN32) || defined(_WIN64))
#include <unistd.h>
#else
#include <io.h>
#endif
#include <time.h>
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"
#include "RmgParallelFft.h"


template  void WriteBGW_VxcEig (int kpt, double *vxc, Kpoint<double> * Kptr);
template  void WriteBGW_VxcEig (int kpt, double *vxc, Kpoint<std::complex<double> > * Kptr);
template <typename KpointType>
void WriteBGW_VxcEig (int kpt, double *vxc, Kpoint<KpointType> * kptr)
{
    KpointType *psi;
    FILE *fhand;

    int nspin = ct.spin_flag +1;
    int ntot_states = ct.num_states * nspin;
    double *vxc_diag = new double[ntot_states];
    int diag_min =ct.vxc_diag_nmin-1, diag_max = ct.vxc_diag_nmax-1;
    int noffdiag = 0;
    int ispin, istate, idx;

    if(diag_max <= diag_min) diag_max = ct.num_states;

    int ndiag = diag_max - diag_min;
    int P0_BASIS = get_P0_BASIS();
    double *vxc_psi = new double[P0_BASIS];


    GetVtotPsi(vxc_psi, vxc, Rmg_G->default_FG_RATIO);
    for(int i = 0; i < ntot_states; i++) vxc_diag[i] = 0.0;

    idx = pct.spinpe * ct.num_states;

    for(istate = diag_min; istate < diag_max; istate++)
    {
        psi = kptr->Kstates[istate].psi;

        for(int i = 0; i < P0_BASIS; i++)
            vxc_diag[idx + istate] += vxc_psi[i] * std::norm(psi[i]);
    }


    MPI_Allreduce(MPI_IN_PLACE, vxc_diag, ntot_states, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    for(int i = 0; i <ntot_states; i++) vxc_diag[i] *= get_vel() *Ha_eV;

    if(pct.gridpe == 0)
    {
        std::string filename("vxc.dat_kpt");
        filename = filename + std::to_string(kpt);

        fhand = fopen((char *)filename.c_str(), "w");
        fprintf(fhand, "%13.9f %13.9f %13.9f %d %d\n", ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ndiag, noffdiag);
        for(ispin = 0; ispin <nspin; ispin++)
            for(istate = diag_min; istate < diag_max; istate++)
            {
                idx = ispin * ct.num_states + istate;
                fprintf(fhand, "%8d %8d %15.9f %15.9f\n", ispin+1, istate+1, vxc_diag[idx],0.0);
            }

        fclose(fhand);
    }
}
