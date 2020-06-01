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


#include <complex>
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Kpoint.h"
#include "transition.h"
#include "Dos.h"

template void OutputDos<double> (Kpoint<double> **);
template void OutputDos<std::complex<double> >(Kpoint<std::complex<double> > **);

template <typename KpointType>
void OutputDos (Kpoint<KpointType> **Kptr)
{

    Kpoint<KpointType> *kptr;
    int nk_tot = ct.num_kpts;
    std::vector<double> eigs;
    eigs.assign(nk_tot * ct.num_states, 0.0);


    for (int ik = 0; ik < ct.num_kpts_pe; ik++)
    {
        kptr = Kptr[ik];
        int ik_glob = pct.kstart + ik;


        for (int is = 0; is < ct.num_states; is++)
        {
            eigs[ik_glob * ct.num_states + is] = kptr->Kstates[is].eig[0] * Ha_eV;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, eigs.data(), nk_tot * ct.num_states, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);

    Dos *Dos_calc = new Dos(ct.kpoint_mesh, ct.kpoint_is_shift, Rmg_L, ct.gaus_broad);
    double Ef_ev = ct.efermi * Ha_eV;
    Dos_calc->tot_dos(nk_tot, ct.num_states, eigs, Ef_ev);

}

