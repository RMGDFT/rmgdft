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
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h> 
#include <unistd.h>
#include <unordered_map>
#include <csignal>
#include <sstream>
#include <iomanip>

#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "transition.h"
#include "Kpoint.h"
#include "InputKey.h"
#include "blas.h"
#include "RmgThread.h"
#include "rmgthreads.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "prototypes_tddft.h"
#include "Exxbase.h"
#include "Neb.h"
#include "Wannier.h"


template void STM_calc (Kpoint<double> **Kptr, double *rho, double bias);
template void STM_calc (Kpoint<std::complex<double>> **Kptr, double *rho, double bias);
template <typename OrbitalType> void STM_calc (Kpoint<OrbitalType> **Kptr, double *rho, double bias)
{

    std::ostringstream streamObj;
    // Set Fixed -Point Notation
    streamObj << std::fixed;
    streamObj << std::setprecision(2);
    //Add double to stream
    streamObj << std::abs(bias);

    std::string bias_string = streamObj.str();
    if(bias < 0.0) bias_string = "m" + bias_string; 
    std::string filename = "stm_bias_" + bias_string + "spin"+std::to_string(pct.spinpe) + ".cube";
    // change the occupation then calculate the charge density
    double Emin = std::min(ct.efermi*Ha_eV, ct.efermi*Ha_eV + bias);
    double Emax = std::max(ct.efermi*Ha_eV, ct.efermi*Ha_eV + bias);
    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++)
    {
        for(int st = 0; st < ct.num_states; st++)
        {

            double eig = Kptr[kpt]->Kstates[st].eig[0]*Ha_eV  ;
            if( eig < Emin)
            {
                Kptr[kpt]->Kstates[st].occupation[0] = 
                    std::exp( - (eig-Emin) * (eig-Emin)/ct.gaus_broad/ct.gaus_broad);
            }
            else if(eig > Emax)
            {
                Kptr[kpt]->Kstates[st].occupation[0] = 
                    std::exp( - (eig-Emax) * (eig-Emin)/ct.gaus_broad/ct.gaus_broad);
            }
            else
            {
                Kptr[kpt]->Kstates[st].occupation[0] = 1.0;
            }

            if(pct.gridpe == 0)printf("\n occ %d %f %f %f", st, ct.efermi * Ha_eV, Kptr[kpt]->Kstates[st].eig[0] * Ha_eV, Kptr[kpt]->Kstates[st].occupation[0]);
        }
    }

    int factor = ct.noncoll_factor * ct.noncoll_factor;

    int ratio = Rmg_G->default_FG_RATIO;
    int FP0_BASIS = Rmg_G->get_P0_BASIS(ratio);

    GetNewRhoPost(Kptr, rho);

    if(!ct.norm_conserving_pp) {
        double *augrho = new double[FP0_BASIS*factor]();
        GetAugRho(Kptr, augrho);
        for(int idx = 0;idx < FP0_BASIS*factor;idx++) rho[idx] += augrho[idx];
        delete [] augrho;
    }

    if (ct.nspin == 2)
        get_rho_oppo (rho,  &rho[FP0_BASIS]);
    if(ct.AFM)
    {
        Rmg_Symm->symmetrize_rho_AFM(rho, &rho[FP0_BASIS]);
    }
    else
    {
        if(Rmg_Symm) Rmg_Symm->symmetrize_grid_object(rho);
        if(ct.noncoll && Rmg_Symm)
            Rmg_Symm->symmetrize_grid_vector(&rho[FP0_BASIS]);
    }

    OutputCubeFile(rho, Rmg_G->default_FG_RATIO, filename);
}

