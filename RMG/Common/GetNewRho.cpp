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
    int max_product = (ct.max_nl + 1) * ct.max_nl / 2;
    double *qtpr;

    if(Verify ("freeze_occupied", true, Kpts[0]->ControlMap)) return;

    double *work = new double[pbasis];
    double *product = new double[max_product];
    OrbitalType *sint = new OrbitalType[2 * ct.max_nl];


    for(int idx = 0;idx < pbasis;idx++)
        work[idx] = 0.0;

    for (int kpt = 0; kpt < ct.num_kpts; kpt++)
    {

        /* Loop over states and accumulate charge */
        for (int istate = 0; istate < nstates; istate++)
        {

            double scale = Kpts[kpt]->Kstates[istate].occupation[0] * ct.kp[kpt].kweight;

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


    if(!ct.norm_conserving_pp) {

        for (int ion = 0; ion < pct.num_nonloc_ions; ion++)
        {
            int gion = pct.nonloc_ions_list[ion];
            
            if (pct.Qidxptrlen[gion])
            {
                
                ION *iptr = &ct.ions[gion];
           
                int nh = ct.sp[iptr->species].nh;
                
                int *ivec = pct.Qindex[gion];
                int ncount = pct.Qidxptrlen[gion];
                double *qnmI = pct.augfunc[gion];

                for (int i=0; i < max_product; i++)
                    product[i] = 0.0;

                for (int kpt = 0; kpt < ct.num_kpts; kpt++)
                {

                    //STATE *sp = ct.kp[kpt].kstate;
                    /* Loop over states and accumulate charge */
                    for (int istate = 0; istate < ct.num_states; istate++)
                    {
                        double t1 = Kpts[kpt]->Kstates[istate].occupation[0] * ct.kp[kpt].kweight;

                        for (int i = 0; i < ct.max_nl; i++)
                        {
                            sint[i] = Kpts[kpt]->newsint_local[ion * ct.num_states * ct.max_nl + istate * ct.max_nl + i];
    //                        sintR[i] =
    //                            pct.newsintR_local[kpt * pct.num_nonloc_ions * ct.num_states * ct.max_nl 
    //                            + ion * ct.num_states * ct.max_nl + istate * ct.max_nl + i];

    //                        if(!ct.is_gamma) {
    //                            sintI[i] =
    //                                pct.newsintI_local[kpt * pct.num_nonloc_ions * ct.num_states * ct.max_nl 
    //                                + ion * ct.num_states * ct.max_nl + istate * ct.max_nl + i];
    //                        }

                        }               /*end for i */

                        int idx = 0;
                        for (int i = 0; i < nh; i++)
                        {
                            for (int j = i; j < nh; j++)
                            {

                                if(i == j) {

                                        product[idx] += t1 * (std::real(sint[i]) * std::real(sint[j]) + std::imag(sint[i]) * std::imag(sint[j]));

                                }
                                else {

                                        product[idx] += 2.0 * t1 * (std::real(sint[i]) * std::real(sint[j]) + std::imag(sint[i]) * std::imag(sint[j]));

                                }
                                idx++;
                            }           /*end for j */
                        }               /*end for i */
                    }                   /*end for istate */
                }                       /*end for kpt */


                int idx = 0;
                for (int i = 0; i < nh; i++)
                {
                    for (int j = i; j < nh; j++)
                    {
                        qtpr = qnmI + idx * ncount;
                        for (int icount = 0; icount < ncount; icount++)
                        {
                            rho[ivec[icount]] += qtpr[icount] * product[idx];
                        }           /*end for icount */
                        idx++;
                    }               /*end for j */
                }                   /*end for i */


            }                       /*end if */

        }                           /*end for ion */

    }

    symmetrize_rho (rho);


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



    delete [] sint;
    delete [] product;
    delete [] work;
}
