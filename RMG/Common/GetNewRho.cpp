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
#include "RmgParallelFft.h"
#include "GlobalSums.h"
#include "Prolong.h"


template void GetNewRho<double>(Kpoint<double> **, double *);
template void GetNewRho<std::complex<double> >(Kpoint<std::complex<double>> **, double *);

template <typename OrbitalType> void GetNewRho(Kpoint<OrbitalType> **Kpts, double *rho)
{

    if(Verify ("freeze_occupied", true, Kpts[0]->ControlMap)) return;

    int pbasis = Kpts[0]->pbasis;
    int nstates = Kpts[0]->nstates;
    int ratio = Rmg_G->default_FG_RATIO;
    int FP0_BASIS = Rmg_G->get_P0_BASIS(ratio);
    Prolong P(ratio, 10, *Rmg_T);

    int dimx = Rmg_G->get_PX0_GRID(ratio);
    int dimy = Rmg_G->get_PY0_GRID(ratio);
    int dimz = Rmg_G->get_PZ0_GRID(ratio);
    int half_dimx = Rmg_G->get_PX0_GRID(1);
    int half_dimy = Rmg_G->get_PY0_GRID(1);
    int half_dimz = Rmg_G->get_PZ0_GRID(1);

    int factor = ct.noncoll_factor * ct.noncoll_factor;
    double *work = new double[FP0_BASIS * factor]();

    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
//#pragma omp parallel
        {
            std::complex<double> psiud;
            double *tarr = new double[FP0_BASIS * factor]();
            OrbitalType *psi_f = new OrbitalType[FP0_BASIS];

            /* Loop over states and accumulate charge */
//#pragma omp for schedule(dynamic)
            for (int istate = 0; istate < nstates; istate++)
            {

                double scale = Kpts[kpt]->Kstates[istate].occupation[0] * Kpts[kpt]->kp.kweight;
                OrbitalType *psi = Kpts[kpt]->Kstates[istate].psi;

//#pragma omp critical(GNR_part1)
                P.prolong(psi_f, psi, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);

                for (int idx = 0; idx < FP0_BASIS; idx++)
                {
                    tarr[idx] += scale * std::norm(psi_f[idx]);
                    if(ct.noncoll)
                    {
                        psiud = 2.0 * psi_f[idx] * std::conj(psi_f[idx + FP0_BASIS]);
                        tarr[idx + 1 * FP0_BASIS] += scale * std::real(psiud);
                        tarr[idx + 2 * FP0_BASIS] += scale * std::imag(psiud);
                        tarr[idx + 3 * FP0_BASIS] += scale * std::norm(psi_f[idx + FP0_BASIS]);
                    }
                }                   /* end for */

            }                       /*end for istate */

//#pragma omp critical(GNR_part2)
            for(int idx = 0; idx < FP0_BASIS * factor; idx++) work[idx] += tarr[idx];

            delete [] psi_f;
            delete [] tarr;
        }
    }                           /*end for kpt */

    MPI_Allreduce(MPI_IN_PLACE, (double *)work, FP0_BASIS, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
    if(ct.noncoll)
    {
        double rho_up, rho_down;
        for(int idx = 0; idx < FP0_BASIS; idx++)
        {
            rho_up = work[idx];
            rho_down = work[idx+3 * FP0_BASIS];
            work[idx] = rho_up + rho_down;
            work[idx+3*FP0_BASIS] = rho_up - rho_down;
        }
    }

    for(int idx = 0; idx < FP0_BASIS * factor; idx++) rho[idx] = work[idx];

#if 0
    /* Interpolate onto fine grid, result will be stored in rho*/
    for(int is = 0; is < factor; is++)
    {
        switch (ct.interp_flag)
        {
            case CUBIC_POLYNOMIAL_INTERPOLATION:
                pack_rho_ctof (&work[is*pbasis], &rho[is*FP0_BASIS]);
                break;
            case PROLONG_INTERPOLATION:
                mg_prolong_MAX10 (&rho[is*FP0_BASIS], &work[is*pbasis], get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);
                break;
            case FFT_INTERPOLATION:
                FftInterpolation (*Kpts[0]->G, &work[is*pbasis], &rho[is*FP0_BASIS], Rmg_G->default_FG_RATIO, ct.sqrt_interpolation);
                break;
            default:

                //Dprintf ("charge interpolation is set to %d", ct.interp_flag);
                rmg_error_handler (__FILE__, __LINE__, "ct.interp_flag is set to an invalid value.");


        }
    }
#endif

    if(!ct.norm_conserving_pp) {
        if(ct.noncoll)
        {
            rmg_error_handler (__FILE__, __LINE__, "ultrasoft with noncollinear is in progress");
        }
        double *augrho = new double[FP0_BASIS]();
        GetAugRho(Kpts, augrho);
        for(int idx = 0;idx < FP0_BASIS;idx++) rho[idx] += augrho[idx];
        delete [] augrho;
    }

    for(int is = 0; is < factor; is++)
        symmetrize_rho (&rho[is*FP0_BASIS]);


    /* Check total charge. */
    ct.tcharge = ZERO;
    for (int idx = 0; idx < FP0_BASIS; idx++)
        ct.tcharge += rho[idx];

    /* ct.tcharge = real_sum_all (ct.tcharge); */
    ct.tcharge = real_sum_all (ct.tcharge, pct.grid_comm);
    ct.tcharge = real_sum_all (ct.tcharge, pct.spin_comm);
    ct.tcharge = ct.tcharge * get_vel_f();

    /* Renormalize charge, there could be some discrepancy because of interpolation */
    double t1 = ct.nel / ct.tcharge;
    for(int i = 0;i < FP0_BASIS * factor;i++) rho[i] *= t1;

    /*Write out normalization constant if needed*/
    double difference = fabs(t1 - 1.0);
    if ((ct.verbose == 1) || (difference > 0.01))
        rmg_printf ("Charge normalization constant: %f\n", t1);


    delete [] work;
}
