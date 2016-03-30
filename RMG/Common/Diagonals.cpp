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


#include "transition.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Solvers.h"
#include "../Headers/prototypes.h"


// Computes the diagonal component of the non-local and S operators
template void Diagonals<double>(Kpoint<double> *);
template void Diagonals<std::complex<double> >(Kpoint<std::complex<double>> *);


template <typename KpointType>
void Diagonals(Kpoint<KpointType> *kptr)
{

    int pbasis = kptr->pbasis;
    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);
    double vel = kptr->L->get_omega() /
                 ((double)(kptr->G->get_NX_GRID(1) * kptr->G->get_NY_GRID(1) * kptr->G->get_NZ_GRID(1)));

    double *dnmI;
    double *qqq;

    for(int idx = 0;idx < pbasis;idx++) kptr->s_diag[idx] = ONE_t;
    for(int idx = 0;idx < pbasis;idx++) kptr->vnl_diag[idx] = ZERO_t;

    if(pct.num_tot_proj == 0)
        return;


    for (int ion = 0; ion < pct.num_nonloc_ions; ion++)
    {

        /*Actual index of the ion under consideration*/
        int proj_index = ion * ct.max_nl;
        int gion = pct.nonloc_ions_list[ion];
        ION *iptr = &ct.ions[gion];
        SPECIES *sp = &ct.sp[iptr->species];

        int nh = sp->nh;

        dnmI = pct.dnmI[gion];
        qqq = pct.qqq[gion];

        for (int i = 0; i < nh; i++)
        {
            int inh = i * nh;
            for (int j = 0; j < nh; j++)
            {
                for(int idx = 0;idx < pbasis;idx++) {
                    kptr->vnl_diag[idx] += dnmI[inh+j] * kptr->nl_weight[(proj_index + i)*pbasis + idx] * kptr->nl_weight[(proj_index + j)*pbasis + idx];
                    kptr->s_diag[idx] += qqq[inh+j] * kptr->nl_weight[(proj_index + i)*pbasis + idx] * kptr->nl_weight[(proj_index + j)*pbasis + idx];
                }

            }
        }
    }


}

