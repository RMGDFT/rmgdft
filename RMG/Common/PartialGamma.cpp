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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

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


/*Set this to 1 to have forces written out part by part*/
/* If you want this , you should also make sure that VERBOSE flag is enabled in
 * nlforce.c*/

template void PartialGamma<double> ( int ion, double * par_gammaR, double * par_omegaR, int nion, int nh,
                     Kpoint<double> **Kptr, int state_start, int state_end, double * sint_derx, double * sint_dery, double * sint_derz );
template void PartialGamma<std::complex<double>> ( int ion, double * par_gammaR, double * par_omegaR, int nion, int nh,
                    Kpoint<std::complex<double>> **Kptr, int state_start, int state_end, std::complex<double> *sint_derx,
                    std::complex<double> *sint_dery, std::complex<double> *sint_derz);


template <typename OrbitalType> void PartialGamma (
        int ion, double * par_gammaR, double * par_omegaR, int nion, int nh,
        Kpoint<OrbitalType> **Kptr, int state_start, int state_end, OrbitalType *sint_derx, OrbitalType *sint_dery, OrbitalType *sint_derz)
{
    int i, j, idx, kidx, istate, size;
    double *gamma_x, *gamma_y, *gamma_z;
    double *omega_x, *omega_y, *omega_z;
    double t1;
    OrbitalType betaxpsiN, betaxpsiM;
    OrbitalType betaxpsiN_der, betaxpsiM_der;
    double beta_pbeta;


    size = (nh * (nh + 1)) / 2;
    gamma_x = par_gammaR;
    gamma_y = gamma_x + size;
    gamma_z = gamma_y + size;

    omega_x = par_omegaR;
    omega_y = omega_x + size;
    omega_z = omega_y + size;


    for (kidx = 0; kidx < ct.num_kpts_pe; kidx++)
    {
        for (istate = state_start; istate < state_end; istate++)
        {
            int istate_local = istate - state_start;
            t1 = Kptr[kidx]->Kstates[istate].occupation[0] * Kptr[kidx]->kweight;



            idx = 0;
            for (i = 0; i < nh; i++)
            {
                for (j = i; j < nh; j++)
                {
                    betaxpsiN = Kptr[kidx]->newsint_local[ istate * pct.num_nonloc_ions * ct.max_nl +
                        nion * ct.max_nl + i];
                    betaxpsiM = Kptr[kidx]->newsint_local[ istate * pct.num_nonloc_ions * ct.max_nl +
                        nion * ct.max_nl + j];

                    betaxpsiN_der = sint_derx[ istate_local * pct.num_nonloc_ions * ct.max_nl +
                        nion * ct.max_nl + i];
                    betaxpsiM_der = sint_derx[ istate_local * pct.num_nonloc_ions * ct.max_nl +
                        nion * ct.max_nl + j];


                    beta_pbeta = std::real( betaxpsiN_der * std::conj(betaxpsiM) +
                            betaxpsiM_der * std::conj(betaxpsiN));
                    gamma_x[idx] += t1 * beta_pbeta;
                    omega_x[idx] += t1 * Kptr[kidx]->Kstates[istate].eig[0] * beta_pbeta;

                    betaxpsiN_der = sint_dery[ istate_local * pct.num_nonloc_ions * ct.max_nl +
                        nion * ct.max_nl + i];
                    betaxpsiM_der = sint_dery[ istate_local * pct.num_nonloc_ions * ct.max_nl +
                        nion * ct.max_nl + j];


                    beta_pbeta = std::real( betaxpsiN_der * std::conj(betaxpsiM) +
                            betaxpsiM_der * std::conj(betaxpsiN));
                    gamma_y[idx] += t1 * beta_pbeta;
                    omega_y[idx] += t1 * Kptr[kidx]->Kstates[istate].eig[0] * beta_pbeta;

                    betaxpsiN_der = sint_derz[ istate_local * pct.num_nonloc_ions * ct.max_nl +
                        nion * ct.max_nl + i];
                    betaxpsiM_der = sint_derz[ istate_local * pct.num_nonloc_ions * ct.max_nl +
                        nion * ct.max_nl + j];


                    beta_pbeta = std::real( betaxpsiN_der * std::conj(betaxpsiM) +
                            betaxpsiM_der * std::conj(betaxpsiN));
                    gamma_z[idx] += t1 * beta_pbeta;
                    omega_z[idx] += t1 * Kptr[kidx]->Kstates[istate].eig[0] * beta_pbeta;


                    ++idx;
                }               /* end for j */
            }                   /*end for i */
        }                       /*end for istate */
    }                           /*end for kidx */

}
