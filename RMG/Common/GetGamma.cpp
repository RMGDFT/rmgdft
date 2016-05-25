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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
#include "GlobalSums.h"


template void GetGamma<double> (double * gammaR, int ion, int nh, Kpoint<double> **Kptr);
template void GetGamma<std::complex<double>> (double * gammaR, int ion, int nh, Kpoint<std::complex<double>> **Kptr);


template <typename OrbitalType> void GetGamma (double * gammaR, int ion, int nh , Kpoint<OrbitalType> **Kptr)
{
    int i, j, idx, kidx, istate;
    double t1;
    OrbitalType sintN, sintM;


    idx = 0;
    for (i = 0; i < nh; i++)
    {
        for (j = i; j < nh; j++)
        {
            gammaR[idx] = 0.0;
            for (kidx = 0; kidx < ct.num_kpts_pe; kidx++)
            {
                for (istate = 0; istate < ct.num_states; istate++)
                {
                    t1 = Kptr[kidx]->Kstates[istate].occupation[0] * Kptr[kidx]->kweight;
                    sintN = Kptr[kidx]->newsint_local[istate *pct.num_nonloc_ions * ct.max_nl +
                        ion * ct.max_nl + i];
                    sintM = Kptr[kidx]->newsint_local[istate *pct.num_nonloc_ions * ct.max_nl +
                        ion * ct.max_nl + j];

                    //gammaR[idx] += t1 * (std::real(sintN) * std::real(sintM) + std::imag(sintN) * std::imag(sintM));
                    gammaR[idx] += t1 * std::real(sintN * std::conj(sintM));
                }               /*end for istate */
            }                   /*end for kidx */
            ++idx;
        }                       /*end for j */
    }                           /*end for i */


    idx = (nh*(nh+1))/2;

    GlobalSums(gammaR, idx, pct.kpsub_comm);

}
