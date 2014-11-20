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


/*Set this to 1 to have forces written out part by part*/
/* If you want this , you should also make sure that VERBOSE flag is enabled in
 * nlforce.c*/

template void PartialGamma<double> ( int ion, double * par_gammaR, double * par_omegaR, int nion, int nh,
                    double * newsintR_x, double * newsintR_y, double * newsintR_z,
                    double * newsintI_x, double * newsintI_y, double * newsintI_z,
                     Kpoint<double> **Kptr);
template void PartialGamma<std::complex<double>> ( int ion, double * par_gammaR, double * par_omegaR, int nion, int nh,
                    double * newsintR_x, double * newsintR_y, double * newsintR_z,
                    double * newsintI_x, double * newsintI_y, double * newsintI_z,
                    Kpoint<std::complex<double>> **Kptr);


template <typename OrbitalType> void PartialGamma (
                    int ion, double * par_gammaR, double * par_omegaR, int nion, int nh,
                    double * newsintR_x, double * newsintR_y, double * newsintR_z,
                    double * newsintI_x, double * newsintI_y, double * newsintI_z,
                    Kpoint<OrbitalType> **Kptr)
{
    int i, j, idx, kidx, istate, size, index;
    double *gamma_x, *gamma_y, *gamma_z;
    double *omega_x, *omega_y, *omega_z;
    double t1;
    OrbitalType betaxpsiN, betaxpsiM;
    double *psixbetaI_x, *psixbetaI_y, *psixbetaI_z;
    double *psixbetaR_x, *psixbetaR_y, *psixbetaR_z;


    size = nh * (nh + 1) / 2;
    gamma_x = par_gammaR;
    gamma_y = gamma_x + size;
    gamma_z = gamma_y + size;

    omega_x = par_omegaR;
    omega_y = omega_x + size;
    omega_z = omega_y + size;

    for (idx = 0; idx < size; idx++)
    {
        gamma_x[idx] = 0.0;
        gamma_y[idx] = 0.0;
        gamma_z[idx] = 0.0;
        omega_x[idx] = 0.0;
        omega_y[idx] = 0.0;
        omega_z[idx] = 0.0;
    }


    for (kidx = 0; kidx < ct.num_kpts; kidx++)
    {
        for (istate = 0; istate < ct.num_states; istate++)
        {
            t1 = Kptr[kidx]->Kstates[istate].occupation[0] * Kptr[kidx]->kweight;


            index =
                ion * ct.num_kpts * ct.num_states * ct.max_nl + kidx * ct.num_states * ct.max_nl +
                istate * ct.max_nl;

            psixbetaR_x = &newsintR_x[index];
            psixbetaR_y = &newsintR_y[index];
            psixbetaR_z = &newsintR_z[index];

            psixbetaI_x = &newsintI_x[index];
            psixbetaI_y = &newsintI_y[index];
            psixbetaI_z = &newsintI_z[index];


            idx = 0;
            for (i = 0; i < nh; i++)
            {
                for (j = i; j < nh; j++)
                {
                    betaxpsiN = Kptr[kidx]->newsint_local[ nion * ct.num_states * ct.max_nl +
                        istate * ct.max_nl + i];
                    betaxpsiM = Kptr[kidx]->newsint_local[ nion * ct.num_states * ct.max_nl +
                        istate * ct.max_nl + j];



                    gamma_x[idx] += t1 * (psixbetaR_x[i] * std::real(betaxpsiM) + psixbetaI_x[i] * std::imag(betaxpsiM)
                                + std::real(betaxpsiN) * psixbetaR_x[j] +
                                std::imag(betaxpsiN) * psixbetaI_x[j]);
                    gamma_y[idx] +=
                        t1 * (psixbetaR_y[i] * std::real(betaxpsiM) + psixbetaI_y[i] * std::imag(betaxpsiM) +
                                std::real(betaxpsiN) * psixbetaR_y[j] + std::imag(betaxpsiN) * psixbetaI_y[j]);
                    gamma_z[idx] +=
                        t1 * (psixbetaR_z[i] * std::real(betaxpsiM) + psixbetaI_z[i] * std::imag(betaxpsiM) +
                                std::real(betaxpsiN) * psixbetaR_z[j] + std::imag(betaxpsiN) * psixbetaI_z[j]);

                    omega_x[idx] +=
                        t1 * Kptr[kidx]->Kstates[istate].eig[0] * (psixbetaR_x[i] * std::real(betaxpsiM) + psixbetaI_x[i] *
                                std::imag(betaxpsiM) +
                                std::real(betaxpsiN) * psixbetaR_x[j] + std::imag(betaxpsiN) * psixbetaI_x[j]);
                    omega_y[idx] +=
                        t1 * Kptr[kidx]->Kstates[istate].eig[0] * (psixbetaR_y[i] * std::real(betaxpsiM) + psixbetaI_y[i] *
                                std::imag(betaxpsiM) +
                                std::real(betaxpsiN) * psixbetaR_y[j] + std::imag(betaxpsiN) * psixbetaI_y[j]);
                    omega_z[idx] +=
                        t1 * Kptr[kidx]->Kstates[istate].eig[0] * (psixbetaR_z[i] * std::real(betaxpsiM) + psixbetaI_z[i] *
                                std::imag(betaxpsiM) +
                                std::real(betaxpsiN) * psixbetaR_z[j] + std::imag(betaxpsiN) * psixbetaI_z[j]);
                    ++idx;
                }               /* end for j */
            }                   /*end for i */
        }                       /*end for istate */
    }                           /*end for kidx */

}
