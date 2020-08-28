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

template void OutputEigenvalues<double> (Kpoint<double> **, int, int);
template void OutputEigenvalues<std::complex<double> >(Kpoint<std::complex<double> > **, int, int);

template <typename KpointType>
void OutputEigenvalues (Kpoint<KpointType> **Kptr, int ikbs, int iscf)
{
    int nk, nspin = (ct.spin_flag + 1);

    int bs = Verify ("calculation_mode", "Band Structure Only", Kptr[0]->ControlMap);

    Kpoint<KpointType> *kptr;
    nk = (bs) ? 1 : ct.num_kpts_pe;

    if(!bs)
    {
        std::vector<double> eigs_up(ct.num_kpts * ct.num_states);
        std::vector<double> eigs_down(ct.num_kpts * ct.num_states);
        std::vector<double> eigs_up_occ(ct.num_kpts * ct.num_states);
        std::vector<double> eigs_down_occ(ct.num_kpts * ct.num_states);
        eigs_up.assign(ct.num_kpts * ct.num_states, 0.0);
        eigs_down.assign(ct.num_kpts * ct.num_states, 0.0);
        eigs_up_occ.assign(ct.num_kpts * ct.num_states, 0.0);
        eigs_down_occ.assign(ct.num_kpts * ct.num_states, 0.0);

        for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++)
        {
            int kpt1 = pct.kstart + kpt;  // absolute kpoint index
            for(int st = 0;st < ct.num_states;st++)
            {   
                if(pct.gridpe==0)
                {
                    eigs_up[kpt1*ct.num_states + st] = Kptr[kpt]->Kstates[st].eig[0];
                    eigs_up_occ[kpt1*ct.num_states + st] = Kptr[kpt]->Kstates[st].occupation[0];
                    if(ct.nspin == 2)
                    {
                        eigs_down[kpt1*ct.num_states + st] = Kptr[kpt]->Kstates[st].eig[1];
                        eigs_down_occ[kpt1*ct.num_states + st] = Kptr[kpt]->Kstates[st].occupation[1];
                    }
                }
            }
        }
        MPI_Allreduce (MPI_IN_PLACE, eigs_up.data(), ct.num_kpts * ct.num_states, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
        MPI_Allreduce (MPI_IN_PLACE, eigs_up_occ.data(), ct.num_kpts * ct.num_states, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
        if(ct.nspin == 2)
        {
            MPI_Allreduce (MPI_IN_PLACE, eigs_down.data(), ct.num_kpts * ct.num_states, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
            MPI_Allreduce (MPI_IN_PLACE, eigs_down_occ.data(), ct.num_kpts * ct.num_states, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
        }

        for(int kpt = 0;kpt < ct.num_kpts;kpt++)
        {
            rmg_printf ("\n\nKOHN SHAM EIGENVALUES [eV] AT K-POINT [%3d]:   %12.6f  %12.6f  %12.6f\n\n",
                    kpt, ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2]);

            for (int idx = 0; idx < nspin; idx++)
            {
                if ( (nspin == 2) && (idx == 0))	
                    rmg_printf("\n------------- SPIN UP ---------------\n\n");
                else if ( (nspin == 2) && (idx == 1))	
                    rmg_printf("\n------------ SPIN DOWN --------------\n\n"); 
                int il = 0;
                for (int is = 0; is < ct.num_states; is++)
                {
                    if (is % 5 == 0)
                        rmg_printf ("[kpt %3d %3d %3d]", kpt, iscf, il++);

                    if(idx == 0)
                        rmg_printf ("   %8.4f [%5.3f]%s",
                            eigs_up[kpt*ct.num_states+is] * Ha_eV, eigs_up_occ[kpt*ct.num_states+is], ((is % 5 == 4) ? "\n" : ""));
                    else
                        rmg_printf ("   %8.4f [%5.3f]%s",
                            eigs_down[kpt*ct.num_states+is] * Ha_eV, eigs_down_occ[kpt*ct.num_states+is], ((is % 5 == 4) ? "\n" : ""));
                }
                rmg_printf ("\n");
           }
           rmg_printf ("\n\n");
        }
    }
    else
    {

        for (int ik = 0; ik < nk; ik++)
        {
            int jk;
            if (bs)
            {
                kptr = Kptr[0];
                jk = ikbs;
            }
            else
            {
                kptr = Kptr[ik];
                jk = ik;
            }


            rmg_printf ("\n\nKOHN SHAM EIGENVALUES [eV] AT K-POINT [%3d]:   %12.6f  %12.6f  %12.6f\n\n",
                    jk, kptr->kp.kpt[0], kptr->kp.kpt[1], kptr->kp.kpt[2]);

            for (int idx = 0; idx < nspin; idx++)
            {
                if ( (nspin == 2) && (idx == 0))	
                    rmg_printf("\n------------- SPIN UP ---------------\n\n");
                else if ( (nspin == 2) && (idx == 1))	
                    rmg_printf("\n------------ SPIN DOWN --------------\n\n"); 
                int il = 0;
                for (int is = 0; is < ct.num_states; is++)
                {
                    if (is % 5 == 0)
                        rmg_printf ("[kpt %3d %3d %3d]", jk, iscf, il++);

                    rmg_printf ("   %8.4f [%5.3f]%s",
                            kptr->Kstates[is].eig[idx] * Ha_eV, kptr->Kstates[is].occupation[idx], ((is % 5 == 4) ? "\n" : ""));
                }
                rmg_printf ("\n");
           }
           rmg_printf ("\n\n");
        }
    }
}

