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
    int ik, jk, nk, is, il, idx, nspin = (ct.spin_flag + 1);

    int bs = Verify ("calculation_mode", "Band Structure Only", Kptr[0]->ControlMap);

    Kpoint<KpointType> *kptr;
    nk = (bs) ? 1 : ct.num_kpts_pe;

    for (ik = 0; ik < nk; ik++)
    {

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
                jk, kptr->kpt[0], kptr->kpt[1], kptr->kpt[2]);

        for (idx = 0; idx < nspin; idx++)
        {
            if ( (nspin == 2) && (idx == 0))	
                rmg_printf("\n------------- SPIN UP ---------------\n\n");
            else if ( (nspin == 2) && (idx == 1))	
                rmg_printf("\n------------ SPIN DOWN --------------\n\n"); 
            il = 0;
            for (is = 0; is < ct.num_states; is++)
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

