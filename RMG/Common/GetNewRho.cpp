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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include <complex>
#include "../Headers/prototypes.h"
#include "RmgParallelFft.h"
#include "GlobalSums.h"


template void GetNewRho<double>(Kpoint<double> **, double *);
template void GetNewRho<std::complex<double> >(Kpoint<std::complex<double>> **, double *);

template <typename OrbitalType> void GetNewRho(Kpoint<OrbitalType> **Kpts, double *rho)
{

    int pbasis = Kpts[0]->pbasis;
    int nstates = Kpts[0]->nstates;
    int max_product = (ct.max_nl + 1) * ct.max_nl / 2;
    double nspin = (ct.spin_flag + 1.0);

    if(Verify ("freeze_occupied", true, Kpts[0]->ControlMap)) return;

    double *work = new double[pbasis];
    double *product = new double[max_product];
    OrbitalType *sint = new OrbitalType[2 * ct.max_nl];


    for(int idx = 0;idx < pbasis;idx++)
        work[idx] = 0.0;

    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {

        /* Loop over states and accumulate charge */
        for (int istate = 0; istate < nstates; istate++)
        {

            double scale = Kpts[kpt]->Kstates[istate].occupation[0] * Kpts[kpt]->kweight;

            OrbitalType *psi = Kpts[kpt]->Kstates[istate].psi;

            for (int idx = 0; idx < pbasis; idx++)
            {
                work[idx] += scale * std::norm(psi[idx]);
            }                   /* end for */

        }                       /*end for istate */

    }                           /*end for kpt */

    GlobalSums(work, pbasis, pct.kpsub_comm);

    /* Interpolate onto fine grid, result will be stored in rho*/
    switch (ct.interp_flag)
    {
        case CUBIC_POLYNOMIAL_INTERPOLATION:
            pack_rho_ctof (work, rho);
            break;
        case BSPLINE_INTERPOLATION:
            bspline_interp_full (work, rho);
            break;
        case PROLONG_INTERPOLATION:
            mg_prolong_MAX10 (rho, work, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);
            break;
        case FFT_INTERPOLATION:
            FftInterpolation (*Kpts[0]->G, work, rho, Rmg_G->default_FG_RATIO);
            break;
        default:

            //Dprintf ("charge interpolation is set to %d", ct.interp_flag);
            rmg_error_handler (__FILE__, __LINE__, "ct.interp_flag is set to an invalid value.");


    }


    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    if(!ct.norm_conserving_pp) {
        double *augrho = new double[FP0_BASIS]();
        GetAugRho(Kpts, augrho);
        for(int idx = 0;idx < FP0_BASIS;idx++) rho[idx] += augrho[idx];
        delete [] augrho;
    }

    symmetrize_rho (rho);


    /* Check total charge. */
    ct.tcharge = ZERO;
    for (int idx = 0; idx < FP0_BASIS; idx++)
        ct.tcharge += rho[idx];

    /* ct.tcharge = real_sum_all (ct.tcharge); */
    ct.tcharge = real_sum_all (ct.tcharge, pct.grid_comm);
    ct.tcharge = ct.tcharge * get_vel_f();

    /* Renormalize charge, there could be some discrepancy because of interpolation */
    double t1 = ct.nel / ct.tcharge / nspin;
//    for(int i = 0;i < FP0_BASIS;i++) rho[i] *= t1;
    
    /*Write out normalization constant if needed*/
    double difference = fabs(t1 - 1.0);
    if ((ct.verbose == 1) || (difference > 0.01))
	rmg_printf ("Charge normalization constant: %f\n", t1);

    /*Update ct.tcharge, do not really recalculate it, just multiply it by normalization constant */
    ct.tcharge *= t1 * nspin;



    delete [] sint;
    delete [] product;
    delete [] work;
}
