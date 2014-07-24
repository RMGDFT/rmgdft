/*
 *
 * Copyright (c) 2014, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
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


template void GetNewRho<double>(Kpoint<double> **, double *);
template void GetNewRho<std::complex<double> >(Kpoint<std::complex<double>> **, double *);

template <typename OrbitalType> void GetNewRho(Kpoint<OrbitalType> **Kpts, double *rho)
{

    int pbasis = Kpts[0]->pbasis;
    int nstates = Kpts[0]->nstates;

    double *work = new double[pbasis];

    for(int idx = 0;idx < pbasis;idx++)
        work[idx] = 0.0;

    for (int kpt = 0; kpt < ct.num_kpts; kpt++)
    {

        /* Loop over states and accumulate charge */
        for (int istate = 0; istate < nstates; istate++)
        {

            double scale = Kpts[kpt]->kstates[istate].occupation[0] * ct.kp[kpt].kweight;

            OrbitalType *psi = Kpts[kpt]->Kstates[istate].psi;

            for (int idx = 0; idx < pbasis; idx++)
            {
                work[idx] += scale * std::norm(psi[idx]);
            }                   /* end for */

        }                       /*end for istate */

    }                           /*end for kpt */

    /* Interpolate onto fine grid, result will be stored in rho*/

    switch (ct.interp_flag)
    {
        case 0:
            pack_rho_ctof (work, rho);
            break;
        case 1:
            bspline_interp_full (work, rho);
            break;
        case 2:
            mg_prolong_MAX10 (rho, work, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);
            break;

        default:

            //Dprintf ("charge interpolation is set to %d", ct.interp_flag);
            rmg_error_handler (__FILE__, __LINE__, "ct.interp_flag is set to an invalid value.");


    }

    if(!ct.is_gamma) {
        symmetrize_rho (rho);
    }


    /* Check total charge. */
    ct.tcharge = ZERO;
    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    for (int idx = 0; idx < FP0_BASIS; idx++)
        ct.tcharge += rho[idx];

    /* ct.tcharge = real_sum_all (ct.tcharge); */
    ct.tcharge = real_sum_all (ct.tcharge, pct.img_comm);
    ct.tcharge = ct.tcharge * get_vel_f();

    /* Renormalize charge, there could be some discrpancy because of interpolation */
    double t1 = ct.nel / ct.tcharge;
    int incx = 1;
    if (pct.imgpe == 0)
        rmg_printf ("\n get_new_rho: Normalization constant for new charge is %f", t1);
    QMD_dscal (FP0_BASIS, t1, rho, incx);

    /*Update ct.tcharge, do not really recalculate it, just mutltiply it by normalization constant */
    ct.tcharge *= t1;



    delete [] work;
}
